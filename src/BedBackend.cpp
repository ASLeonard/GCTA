/*
   BedBackend — coroutine-based genotype I/O for PLINK BED files.

   Replaces:
     preGenoDouble_bed()   → constructor
     readGeno_bed()        → stream() coroutine body
     endGenoDouble_bed()   → destructor

   Resource safety: masks are owned by this object.  If stream() is
   abandoned early (e.g. an exception from the callback), the coroutine
   frame is destroyed by exec::task's destructor, and then unique_ptr
   destroys the BedBackend, which frees the masks.  No manual
   endGenoDouble() call is ever needed.
*/

#include "GenoBackendFactory.h"
#include "Geno.h"                               // full definition for friend access
#include "Logger.h"
#include "mem.hpp"
#include "third_party/Pgenlib/PgenReader.h"

#include <stdexec/execution.hpp>
#include <exec/task.hpp>
#include <exec/static_thread_pool.hpp>

#include <algorithm>
#include <numeric>
#include <vector>

// --------------------------------------------------------------------------
// Helper: aligned-allocated mask buffer with posix_mem_free deleter.
// --------------------------------------------------------------------------
namespace {
struct PosixFreeDeleter {
    void operator()(uintptr_t *p) const noexcept { posix_mem_free(p); }
};
using MaskBuf = std::unique_ptr<uintptr_t, PosixFreeDeleter>;

MaskBuf allocMask(std::size_t nElements)
{
    uintptr_t *p = nullptr;
    if (posix_memalign(reinterpret_cast<void **>(&p), 32,
                       nElements * sizeof(uintptr_t)) != 0)
        p = nullptr;
    return MaskBuf{p};
}
} // anonymous namespace

// --------------------------------------------------------------------------
// BedBackend
// --------------------------------------------------------------------------
class BedBackend : public GenoBackend {
public:
    explicit BedBackend(Geno &geno);
    ~BedBackend() override;

    exec::task<void> stream(GenoScheduler               io_sched,
                            GenoScheduler               cpu_sched,
                            const std::vector<uint32_t> &extractIndex,
                            Callback                    callback) override;

private:
    Geno &g_;

    // Owned mask buffers (lifetime = this backend).
    // The raw pointers in Geno (g_.keepMaskPtr etc.) point into these.
    MaskBuf keepMask_;
    MaskBuf keepMaskInter_;
    MaskBuf sexMask_;
    MaskBuf sexMaskInter_;
    MaskBuf maleMask_;
    MaskBuf maleMaskInter_;
};

// --------------------------------------------------------------------------
BedBackend::BedBackend(Geno &geno) : g_(geno)
{
    // ── Mirrors preGenoDouble_bed() ──────────────────────────────────────
    g_.hasInfo        = false;

    g_.compressFormats.clear();
    g_.rawCountSamples.clear();
    g_.rawCountSNPs.clear();

    for (std::size_t i = 0; i < g_.geno_files.size(); ++i) {
        MarkerParam p = g_.marker->getMarkerParams(static_cast<int>(i));
        g_.compressFormats.push_back(p.compressFormat);
        g_.rawCountSamples.push_back(g_.rawSampleCT);
        g_.rawCountSNPs.push_back(p.rawCountSNP);
    }

    uint32_t rawCT    = g_.rawSampleCT;
    g_.bedRawGenoBuf1PtrSize = PgenReader::GetGenoBufPtrSize(rawCT);
    g_.maskPtrSize           = static_cast<int>(PgenReader::GetSubsetMaskSize(rawCT));

    std::size_t nMask = static_cast<std::size_t>(g_.maskPtrSize);
    keepMask_       = allocMask(nMask);
    keepMaskInter_  = allocMask(nMask);
    sexMask_        = allocMask(nMask);
    sexMaskInter_   = allocMask(nMask);
    maleMask_       = allocMask(nMask);
    maleMaskInter_  = allocMask(nMask);

    if (!keepMask_ || !keepMaskInter_ ||
        !sexMask_  || !sexMaskInter_  ||
        !maleMask_ || !maleMaskInter_) {
        LOGGER.e(0, "BedBackend: can't allocate mask buffers.");
    }

    // Expose raw pointers so that Geno's getGenoDouble_bed can use them.
    g_.keepMaskPtr      = keepMask_.get();
    g_.keepMaskInterPtr = keepMaskInter_.get();
    g_.sexMaskPtr       = sexMask_.get();
    g_.sexMaskInterPtr  = sexMaskInter_.get();
    g_.maleMaskPtr      = maleMask_.get();
    g_.maleMaskInterPtr = maleMaskInter_.get();

    PgenReader::SetSampleSubsets(g_.sampleKeepIndex, rawCT,
                                 g_.keepMaskPtr, g_.keepMaskInterPtr);
    PgenReader::SetSampleSubsets(g_.keepSexIndex,  rawCT,
                                 g_.sexMaskPtr,  g_.sexMaskInterPtr);
    PgenReader::SetSampleSubsets(g_.keepMaleIndex, rawCT,
                                 g_.maleMaskPtr, g_.maleMaskInterPtr);
}

// --------------------------------------------------------------------------
BedBackend::~BedBackend()
{
    // ── Mirrors endGenoDouble_bed() ──────────────────────────────────────
    // MaskBuf unique_ptrs free the memory; null out Geno's raw pointers
    // so that any accidental post-destruction use fails fast.
    g_.keepMaskPtr      = nullptr;
    g_.keepMaskInterPtr = nullptr;
    g_.sexMaskPtr       = nullptr;
    g_.sexMaskInterPtr  = nullptr;
    g_.maleMaskPtr      = nullptr;
    g_.maleMaskInterPtr = nullptr;
}

// --------------------------------------------------------------------------
// NOTE: callback is taken by value so it lives in the coroutine frame and
// remains valid across every co_await suspension point.  Do NOT change to a
// reference — it would dangle the moment the caller's scope is suspended.
exec::task<void> BedBackend::stream(GenoScheduler               io_sched,
                                    GenoScheduler               cpu_sched,
                                    const std::vector<uint32_t> &extractIndex,
                                    Callback                    callback)
try {
    // ── Transfer to the I/O thread ───────────────────────────────────────
    co_await stdexec::schedule(io_sched);

    // ── Build rawIndices (extract → raw marker index) ────────────────────
    const std::vector<uint32_t> &raw_marker_index = g_.marker->get_extract_index();
    std::vector<uint32_t> rawIndices(extractIndex.size());
    std::transform(extractIndex.begin(), extractIndex.end(),
                   rawIndices.begin(),
                   [&raw_marker_index](std::size_t pos) {
                       return raw_marker_index[pos];
                   });

    const uint32_t numMarker    = static_cast<uint32_t>(extractIndex.size());
    const int      stride       = g_.bedRawGenoBuf1PtrSize;
    const int      markerBlock  = g_.numMarkerBlock;

    uint32_t finishedMarker = 0;
    int      preFileIndex   = -1;
    int      fileIndex      = 0;
    bool     chr_ends       = false;
    uint8_t  isSexXY        = 0;
    PgenReader reader;

    // ── Main I/O loop ────────────────────────────────────────────────────
    while (finishedMarker < numMarker) {
        uint32_t nextSize = g_.marker->getNextSize(
            rawIndices, finishedMarker, markerBlock,
            fileIndex, chr_ends, isSexXY);
        if (nextSize == 0) break;

        // Build a block on the I/O thread.
        GenoBlock block;
        block.numMarkers = nextSize;
        block.isSexXY   = isSexXY;
        block.fileIndex  = fileIndex;
        block.buf.resize(static_cast<std::size_t>(stride) * nextSize);
        block.extractIndex.assign(
            extractIndex.begin() + finishedMarker,
            extractIndex.begin() + finishedMarker + nextSize);

        if (preFileIndex != fileIndex) {
            reader.Load(g_.geno_files[fileIndex],
                        &g_.rawCountSamples[fileIndex],
                        &g_.rawCountSNPs[fileIndex],
                        g_.sampleKeepIndex);
            preFileIndex = fileIndex;
        }

        for (uint32_t i = 0; i < nextSize; ++i) {
            int rawIndex  = static_cast<int>(rawIndices[finishedMarker + i]);
            int lag_index = rawIndex - g_.baseIndexLookup[fileIndex];
            reader.ReadRawFullHard(
                block.buf.data() + static_cast<std::size_t>(i) * stride,
                lag_index);
        }

        // ── Transfer to CPU thread — callback does compute ────────────────
        co_await stdexec::schedule(cpu_sched);
        callback(block);

        // ── Return to I/O thread for next block ───────────────────────────
        co_await stdexec::schedule(io_sched);

        finishedMarker += nextSize;
    }
    // Destructor handles mask cleanup.
} catch (...) {
    // Re-propagate through exec::task's error channel so that
    // stdexec::sync_wait re-throws on the loopDouble call-site.
    // RAII: ~BedBackend() still runs (masks freed) before the exception
    // surfaces, because backend is destroyed at the end of loopDouble().
    throw;
}

// --------------------------------------------------------------------------
// Factory
// --------------------------------------------------------------------------
std::unique_ptr<GenoBackend> makeBedBackend(Geno &geno)
{
    return std::make_unique<BedBackend>(geno);
}
