/*
   VcfBackend — coroutine-based genotype I/O for VCF / BCF files.

   Format support:
     BCF (CSI index)          — preferred, fast random access
     VCF.gz (CSI index)       — supported
     VCF.gz (TBI index)       — supported with a warning (may fail on long chrs)
     VCF (plain, no index)    — rejected at open time

   Implements the same GenoBackend interface as BedBackend / BgenBackend /
   PgenBackend, feeding the same GenoBlock callback used by loopDouble().

   Buffer layout per marker (vcfRawGenoBuf1PtrSize uintptr_t words):
     keepSampleCT bytes (uint8_t), values:
       0, 1, 2  — allele count of the ALT allele (a1 / effect allele)
       255      — missing genotype
   Samples are ordered by kept-extract position (identical to the order in
   which sampleKeepIndex was built by Pheno).

   Processing (AF computation, centering, standardising, missing mask) is
   done by Geno::getGenoDouble_vcf(), exactly mirroring getGenoDouble_bed().

   Seek strategy: one bcf_itr_querys() region per block (all markers in a
   block share the same chromosome — guaranteed by Marker::getNextSize()).
   This decompresses each BGZF block at most once per marker block.
*/

#include "GenoBackendFactory.h"
#include "Geno.h"
#include "Logger.h"

#include <stdexec/execution.hpp>
#include <exec/task.hpp>
#include <exec/static_thread_pool.hpp>

#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <vector>

#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"

// --------------------------------------------------------------------------
// SampleRemap
//
// Maps VCF-column index → extract-order position in the kept-sample output.
// Entry is -1 if that VCF sample is not in the kept set.
// --------------------------------------------------------------------------
struct SampleRemap {
    std::vector<int32_t> vcfToExtract;  ///< length = n_samples in VCF header
    int32_t              nKept = 0;
};

static SampleRemap buildRemap(bcf_hdr_t *hdr, const std::vector<uint32_t> &keepIndex,
                               Pheno *pheno)
{
    SampleRemap remap;
    int nVcfSamples = bcf_hdr_nsamples(hdr);
    remap.vcfToExtract.assign(static_cast<std::size_t>(nVcfSamples), -1);

    // Build lookup: IID → position in keepIndex.
    // Try IID alone first; also register FID_IID as fallback key.
    std::unordered_map<std::string, int32_t> idToKeepPos;
    idToKeepPos.reserve(keepIndex.size() * 2);
    for (int32_t i = 0; i < static_cast<int32_t>(keepIndex.size()); ++i) {
        uint32_t rawIdx = keepIndex[i];
        const std::string &iid = pheno->getRawIID(rawIdx);
        const std::string &fid = pheno->getRawFID(rawIdx);
        // IID lookup (preferred for standard VCF)
        idToKeepPos.emplace(iid, i);
        // FID_IID lookup (for VCFs generated from PLINK)
        idToKeepPos.emplace(fid + "_" + iid, i);
    }

    for (int v = 0; v < nVcfSamples; ++v) {
        const char *vcfID = hdr->samples[v];
        auto it = idToKeepPos.find(vcfID);
        if (it != idToKeepPos.end()) {
            remap.vcfToExtract[static_cast<std::size_t>(v)] = it->second;
            ++remap.nKept;
        }
    }

    if (remap.nKept == 0)
        LOGGER.e(0, "VcfBackend: no samples in the VCF match the phenotype file. "
                    "Check sample ID format (IID vs FID_IID).");
    if (remap.nKept < static_cast<int32_t>(keepIndex.size())) {
        LOGGER << "Warning: VCF is missing "
               << (keepIndex.size() - static_cast<std::size_t>(remap.nKept))
               << " kept sample(s) from the phenotype file." << std::endl;
    }
    return remap;
}

// --------------------------------------------------------------------------
// Normalise one VCF record's GT field into the uint8_t output buffer.
// Returns false if the record should be skipped (no GT, or multiallelic when
// policy == Skip).
// --------------------------------------------------------------------------
enum class MultiAllelicPolicy : uint8_t { Skip, SplitRef };

static bool normaliseGenotypes(bcf_hdr_t *hdr, bcf1_t *rec,
                                const SampleRemap &remap,
                                MultiAllelicPolicy policy,
                                uint8_t *outBuf,   ///< length >= remap.nKept
                                uint32_t keepCT)
{
    bcf_unpack(rec, BCF_UN_FMT);

    int32_t *gt_arr = nullptr;
    int      n_gt   = 0;
    int      n      = bcf_get_genotypes(hdr, rec, &gt_arr, &n_gt);
    if (n <= 0) { free(gt_arr); return false; }

    if (rec->n_allele > 2 && policy == MultiAllelicPolicy::Skip) {
        free(gt_arr);
        return false;
    }

    int nVcf   = bcf_hdr_nsamples(hdr);
    int ploidy = (nVcf > 0) ? (n_gt / nVcf) : 2;

    // Pre-fill output with 255 (missing) so samples not in the VCF stay missing.
    std::fill(outBuf, outBuf + keepCT, static_cast<uint8_t>(255));

    for (int v = 0; v < nVcf; ++v) {
        int32_t extractPos = remap.vcfToExtract[static_cast<std::size_t>(v)];
        if (extractPos < 0) continue;

        int32_t a0 = gt_arr[v * ploidy + 0];
        int32_t a1_gt = (ploidy > 1) ? gt_arr[v * ploidy + 1] : bcf_int32_vector_end;

        if (bcf_gt_is_missing(a0) || a1_gt == bcf_int32_vector_end ||
            (ploidy > 1 && bcf_gt_is_missing(a1_gt))) {
            outBuf[extractPos] = 255;
        } else {
            int al0 = bcf_gt_allele(a0);
            int al1 = (ploidy > 1) ? bcf_gt_allele(a1_gt) : al0;
            int g;
            if (rec->n_allele > 2) {
                // SplitRef: any non-ref allele counts as 1 alt dosage
                g = (al0 > 0 ? 1 : 0) + (al1 > 0 ? 1 : 0);
            } else {
                g = al0 + al1;
            }
            // Clamp to [0,2] for safety
            outBuf[extractPos] = static_cast<uint8_t>(g < 0 ? 0 : (g > 2 ? 2 : g));
        }
    }

    free(gt_arr);
    return true;
}

// --------------------------------------------------------------------------
// VcfBackend
// --------------------------------------------------------------------------
class VcfBackend : public GenoBackend {
public:
    explicit VcfBackend(Geno &geno);
    ~VcfBackend() override;

    exec::task<void> stream(GenoScheduler               io_sched,
                            GenoScheduler               cpu_sched,
                            const std::vector<uint32_t> &extractIndex,
                            Callback                    callback) override;

private:
    Geno &g_;

    struct VcfHandle {
        htsFile   *fp        = nullptr;
        bcf_hdr_t *hdr       = nullptr;
        hts_idx_t *idx       = nullptr;  ///< BCF path: CSI index
        tbx_t     *tbx       = nullptr;  ///< VCF.gz path: TBI or CSI via tabix API
        bool       is_vcf_gz = false;
        SampleRemap remap;
    };
    std::vector<VcfHandle> handles_;

    MultiAllelicPolicy policy_ = MultiAllelicPolicy::Skip;
};

// --------------------------------------------------------------------------
VcfBackend::VcfBackend(Geno &geno) : g_(geno)
{
    g_.hasInfo = false;

    g_.compressFormats.clear();
    g_.rawCountSamples.clear();
    g_.rawCountSNPs.clear();

    // Compute buf stride: keepSampleCT bytes padded up to uintptr_t words.
    const uint32_t keepCT    = g_.keepSampleCT;
    const std::size_t stride =
        (static_cast<std::size_t>(keepCT) + sizeof(uintptr_t) - 1) / sizeof(uintptr_t);
    g_.vcfRawGenoBuf1PtrSize = static_cast<int>(stride);

    handles_.reserve(g_.geno_files.size());
    for (std::size_t i = 0; i < g_.geno_files.size(); ++i) {
        const std::string &path = g_.geno_files[i];

        htsFile *fp = hts_open(path.c_str(), "r");
        if (!fp)
            LOGGER.e(0, "VcfBackend: cannot open [" + path + "]");

        // Check format
        const htsFormat &fmt = fp->format;
        if (fmt.format == vcf && fmt.compression == no_compression) {
            hts_close(fp);
            LOGGER.e(0, "plain uncompressed VCF [" + path + "] is not supported. "
                        "Convert to BCF or bgzipped VCF:\n"
                        "  bcftools view -O b -o out.bcf in.vcf && bcftools index out.bcf");
        }

        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        if (!hdr) { hts_close(fp); LOGGER.e(0, "VcfBackend: cannot read header from [" + path + "]"); }

        bool       is_gz = (fmt.format == vcf && fmt.compression == bgzf);
        hts_idx_t *idx   = nullptr;
        tbx_t     *tbx   = nullptr;
        if (is_gz) {
            // VCF.gz: load via tabix API (supports both TBI and CSI).
            tbx = tbx_index_load2(path.c_str(), nullptr);
            if (!tbx) {
                bcf_hdr_destroy(hdr); hts_close(fp);
                LOGGER.e(0, "VcfBackend: no index for [" + path + "]. "
                            "Create a TBI index: tabix -p vcf " + path
                            + "\n  or a CSI index: bcftools index --csi " + path);
            }
        } else {
            // BCF: load CSI index for bcf_itr_querys / bcf_itr_next.
            idx = bcf_index_load2(path.c_str(), nullptr);
            if (!idx) {
                bcf_hdr_destroy(hdr); hts_close(fp);
                LOGGER.e(0, "VcfBackend: no index for [" + path + "]. "
                            "Create with: bcftools index --csi " + path);
            }
        }

        SampleRemap remap = buildRemap(hdr, g_.sampleKeepIndex, g_.pheno);

        // Populate Geno's format arrays for baseIndexLookup computation.
        MarkerParam p = g_.marker->getMarkerParams(static_cast<int>(i));
        g_.compressFormats.push_back(0);
        g_.rawCountSamples.push_back(static_cast<uint32_t>(bcf_hdr_nsamples(hdr)));
        g_.rawCountSNPs.push_back(p.rawCountSNP);

        handles_.push_back({fp, hdr, idx, tbx, is_gz, std::move(remap)});
    }

    LOGGER.i(0, "VCF/BCF backend: multiallelic sites will be skipped (default).");
}

// --------------------------------------------------------------------------
VcfBackend::~VcfBackend()
{
    for (auto &h : handles_) {
        if (h.is_vcf_gz) {
            if (h.tbx) tbx_destroy(h.tbx);
        } else {
            if (h.idx) hts_idx_destroy(h.idx);
        }
        if (h.hdr) bcf_hdr_destroy(h.hdr);
        if (h.fp)  hts_close(h.fp);
    }
}

// --------------------------------------------------------------------------
exec::task<void> VcfBackend::stream(GenoScheduler               io_sched,
                                     GenoScheduler               cpu_sched,
                                     const std::vector<uint32_t> &extractIndex,
                                     Callback                    callback)
try {
    co_await stdexec::schedule(io_sched);

    const std::vector<uint32_t> &raw_marker_index = g_.marker->get_extract_index();
    std::vector<uint32_t> rawIndices(extractIndex.size());
    std::transform(extractIndex.begin(), extractIndex.end(),
                   rawIndices.begin(),
                   [&raw_marker_index](std::size_t pos) {
                       return raw_marker_index[pos];
                   });

    const uint32_t numMarker   = static_cast<uint32_t>(extractIndex.size());
    const int      stride      = g_.vcfRawGenoBuf1PtrSize;
    const int      markerBlock = g_.numMarkerBlock;
    const uint32_t keepCT      = g_.keepSampleCT;

    uint32_t finishedMarker = 0;
    int      fileIndex      = 0;
    bool     chr_ends       = false;
    uint8_t  isSexXY        = 0;

    bcf1_t *rec = bcf_init();
    if (!rec)
        LOGGER.e(0, "VcfBackend: bcf_init() failed");

    while (finishedMarker < numMarker) {
        uint32_t nextSize = g_.marker->getNextSize(
            rawIndices, finishedMarker, markerBlock,
            fileIndex, chr_ends, isSexXY);
        if (nextSize == 0) break;

        // ── Build GenoBlock ──────────────────────────────────────────────
        GenoBlock block;
        block.numMarkers = nextSize;
        block.isSexXY    = isSexXY;
        block.fileIndex  = fileIndex;
        block.buf.assign(static_cast<std::size_t>(stride) * nextSize, 0);
        block.extractIndex.assign(
            extractIndex.begin() + finishedMarker,
            extractIndex.begin() + finishedMarker + nextSize);

        // ── Build position → slot map for this block ─────────────────────
        // Use multimap to handle duplicate positions (rare but possible).
        // key = 1-based position (matching Marker::pd storage).
        std::unordered_multimap<uint32_t, uint32_t> posToSlot;
        posToSlot.reserve(nextSize);
        uint32_t firstPos = UINT32_MAX, lastPos = 0;
        for (uint32_t i = 0; i < nextSize; ++i) {
            uint32_t rawIdx = rawIndices[finishedMarker + i];
            uint32_t bp     = g_.marker->getRawBp(rawIdx);
            posToSlot.emplace(bp, i);
            if (bp < firstPos) firstPos = bp;
            if (bp > lastPos)  lastPos  = bp;
        }

        // ── Chromosome string for region query ───────────────────────────
        uint32_t firstRawIdx = rawIndices[finishedMarker];
        const std::string &chrStr = g_.marker->getRawChr(firstRawIdx);

        VcfHandle &h = handles_[static_cast<std::size_t>(fileIndex)];

        // ── One region query per block (§8 region-batched approach) ──────
        std::string region = chrStr + ":" + std::to_string(firstPos)
                             + "-" + std::to_string(lastPos);
        if (h.is_vcf_gz) {
            // VCF.gz: tbx_itr_next reads a text line; vcf_parse fills rec.
            kstring_t kstr = KS_INITIALIZE;
            hts_itr_t *itr = tbx_itr_querys(h.tbx, region.c_str());
            if (itr) {
                while (tbx_itr_next(h.fp, h.tbx, itr, &kstr) >= 0) {
                    if (vcf_parse(&kstr, h.hdr, rec) < 0) continue;
                    uint32_t recPos = static_cast<uint32_t>(rec->pos + 1);
                    auto range = posToSlot.equal_range(recPos);
                    for (auto it = range.first; it != range.second; ++it) {
                        uint32_t slot = it->second;
                        uint8_t *dest = reinterpret_cast<uint8_t *>(
                            block.buf.data() + static_cast<std::size_t>(slot) * stride);
                        std::fill(dest, dest + keepCT, static_cast<uint8_t>(255));
                        normaliseGenotypes(h.hdr, rec, h.remap, policy_, dest, keepCT);
                    }
                }
                hts_itr_destroy(itr);
            }
            ks_free(&kstr);
        } else {
            // BCF: bcf_itr_next calls bcf_readrec which reads BCF binary correctly.
            hts_itr_t *itr = bcf_itr_querys(h.idx, h.hdr, region.c_str());
            if (itr) {
                while (bcf_itr_next(h.fp, itr, rec) >= 0) {
                    uint32_t recPos = static_cast<uint32_t>(rec->pos + 1);
                    auto range = posToSlot.equal_range(recPos);
                    for (auto it = range.first; it != range.second; ++it) {
                        uint32_t slot = it->second;
                        uint8_t *dest = reinterpret_cast<uint8_t *>(
                            block.buf.data() + static_cast<std::size_t>(slot) * stride);
                        std::fill(dest, dest + keepCT, static_cast<uint8_t>(255));
                        normaliseGenotypes(h.hdr, rec, h.remap, policy_, dest, keepCT);
                    }
                }
                hts_itr_destroy(itr);
            }
        }

        // ── Hand off to CPU thread for callback ──────────────────────────
        co_await stdexec::schedule(cpu_sched);
        callback(block);
        co_await stdexec::schedule(io_sched);

        finishedMarker += nextSize;
    }

    bcf_destroy(rec);
} catch (...) {
    // Re-propagate so stdexec::sync_wait re-throws on the loopDouble call-site.
    // RAII: ~VcfBackend() closes all handles before the exception surfaces.
    throw;
}

// --------------------------------------------------------------------------
std::unique_ptr<GenoBackend> makeVcfBackend(Geno &geno)
{
    return std::make_unique<VcfBackend>(geno);
}
