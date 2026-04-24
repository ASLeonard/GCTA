/*
   BgenBackend — coroutine-based genotype I/O for Oxford BGEN files.

   Replaces:
     preGenoDouble_bgen()  → constructor
     readGeno_bgen()       → stream() coroutine body
     endGenoDouble_bgen()  → destructor (trivial — no allocated masks)

   BGEN reads are position-based (fseek), so we need open file handles.
   The handles are opened once per loopDouble call (in the constructor)
   and closed in the destructor.  The raw compressed blob for each marker
   is copied into GenoBlock.buf; decompression happens later in
   Geno::getGenoDouble_bgen() exactly as it did before.
*/

#include "GenoBackendFactory.h"
#include "Geno.h"
#include "Logger.h"

#include <stdexec/execution.hpp>
#include <exec/task.hpp>
#include <exec/static_thread_pool.hpp>

#include <algorithm>
#include <cstdio>
#include <vector>

// --------------------------------------------------------------------------
class BgenBackend : public GenoBackend {
public:
    explicit BgenBackend(Geno &geno);
    ~BgenBackend() override;

    exec::task<void> stream(GenoScheduler               io_sched,
                            GenoScheduler               cpu_sched,
                            const std::vector<uint32_t> &extractIndex,
                            Callback                    callback) override;

private:
    Geno                &g_;
    std::vector<FILE *>  fileHandles_; ///< one per geno_file
};

// --------------------------------------------------------------------------
BgenBackend::BgenBackend(Geno &geno) : g_(geno)
{
    // ── Mirrors preGenoDouble_bgen() ─────────────────────────────────────
    g_.hasInfo = true;

    g_.compressFormats.clear();
    g_.rawCountSamples.clear();
    g_.rawCountSNPs.clear();

    for (std::size_t i = 0; i < g_.geno_files.size(); ++i) {
        MarkerParam p = g_.marker->getMarkerParams(static_cast<int>(i));
        if (p.rawCountSample != g_.pheno->count_raw()) {
            LOGGER.e(0, "inconsistent sample sizes between the .bgen file ["
                        + g_.geno_files[i]
                        + "] and the .sample file (specified by --sample).");
        }
        g_.compressFormats.push_back(p.compressFormat);
        g_.rawCountSamples.push_back(p.rawCountSample);
        g_.rawCountSNPs.push_back(p.rawCountSNP);
    }

    g_.bgenRawGenoBuf1PtrSize = g_.marker->getMaxGenoMarkerUptrSize();

    // Open file handles.
    fileHandles_.reserve(g_.geno_files.size());
    for (const auto &path : g_.geno_files) {
        FILE *f = fopen(path.c_str(), "rb");
        if (!f)
            LOGGER.e(0, "BgenBackend: can't open [" + path + "] to read.");
        fileHandles_.push_back(f);
    }
}

// --------------------------------------------------------------------------
BgenBackend::~BgenBackend()
{
    for (FILE *f : fileHandles_)
        if (f) fclose(f);
}

// --------------------------------------------------------------------------
// NOTE: callback is taken by value so it lives in the coroutine frame and
// remains valid across every co_await suspension point.  Do NOT change to a
// reference — it would dangle the moment the caller's scope is suspended.
exec::task<void> BgenBackend::stream(GenoScheduler               io_sched,
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
    const int      stride      = g_.bgenRawGenoBuf1PtrSize;
    const int      markerBlock = g_.numMarkerBlock;

    uint32_t finishedMarker = 0;
    int      fileIndex      = 0;
    bool     chr_ends       = false;
    uint8_t  isSexXY        = 0;

    while (finishedMarker < numMarker) {
        uint32_t nextSize = g_.marker->getNextSize(
            rawIndices, finishedMarker, markerBlock,
            fileIndex, chr_ends, isSexXY);
        if (nextSize == 0) break;

        GenoBlock block;
        block.numMarkers = nextSize;
        block.isSexXY   = isSexXY;
        block.fileIndex  = fileIndex;
        block.buf.resize(static_cast<std::size_t>(stride) * nextSize);
        block.extractIndex.assign(
            extractIndex.begin() + finishedMarker,
            extractIndex.begin() + finishedMarker + nextSize);

        FILE *bgenFile = fileHandles_[fileIndex];
        for (uint32_t i = 0; i < nextSize; ++i) {
            int rawIndex = static_cast<int>(rawIndices[finishedMarker + i]);

            uint64_t pos, size;
            g_.marker->getStartPosSize(rawIndex, pos, size);
            fseek(bgenFile, static_cast<long>(pos), SEEK_SET);

            uintptr_t *dest = block.buf.data()
                              + static_cast<std::size_t>(i) * stride;
            if (fread(dest, sizeof(char), size, bgenFile) != size) {
                int lag_index = rawIndex - g_.baseIndexLookup[fileIndex];
                LOGGER.e(0, "can't read " + std::to_string(lag_index)
                            + "th SNP in [" + g_.geno_files[fileIndex] + "].");
            }
        }

        co_await stdexec::schedule(cpu_sched);
        callback(block);
        co_await stdexec::schedule(io_sched);

        finishedMarker += nextSize;
    }
} catch (...) {
    // Re-propagate through exec::task's error channel so that
    // stdexec::sync_wait re-throws on the loopDouble call-site.
    // RAII: ~BgenBackend() closes all file handles before the exception
    // surfaces.
    throw;
}

// --------------------------------------------------------------------------
std::unique_ptr<GenoBackend> makeBgenBackend(Geno &geno)
{
    return std::make_unique<BgenBackend>(geno);
}
