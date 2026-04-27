/*
   PgenBackend — coroutine-based genotype I/O for PLINK PGEN files.

   This is the PROPER implementation.  Previously PGEN was silently
   aliased to the BED code path (preGenoDouble_bed / readGeno_bed /
   getGenoDouble_bed), which meant the PGEN-specific dosage buffer
   layout was never set up and readGeno_pgen / getGenoDouble_pgen were
   dead code.

   Replaces:
     preGenoDouble_pgen()  → constructor
     readGeno_pgen()       → stream() coroutine body (uses ReadDosage)
     endGenoDouble_pgen()  → destructor (trivial)

   Processing (getGenoDouble_pgen) is unchanged and remains in Geno.cpp.
*/

#include "GenoBackendFactory.h"
#include "Geno.h"
#include "Logger.h"
#include "mem.hpp"
#include "third_party/Pgenlib/PgenReader.h"

#include <stdexec/execution.hpp>
#include <exec/task.hpp>
#include <exec/static_thread_pool.hpp>

#include <algorithm>
#include <vector>

// --------------------------------------------------------------------------
class PgenBackend : public GenoBackend {
public:
    explicit PgenBackend(Geno &geno);
    ~PgenBackend() override = default;

    exec::task<void> stream(GenoScheduler               io_sched,
                            GenoScheduler               cpu_sched,
                            const std::vector<uint32_t> &extractIndex,
                            Callback                    callback) override;

private:
    Geno &g_;
};

// --------------------------------------------------------------------------
PgenBackend::PgenBackend(Geno &geno) : g_(geno)
{
    // ── Mirrors preGenoDouble_pgen() ─────────────────────────────────────
    g_.hasInfo = false;

    g_.compressFormats.clear();
    g_.rawCountSamples.clear();
    g_.rawCountSNPs.clear();

    for (std::size_t i = 0; i < g_.geno_files.size(); ++i) {
        MarkerParam p = g_.marker->getMarkerParams(static_cast<int>(i));
        g_.compressFormats.push_back(p.compressFormat);
        g_.rawCountSamples.push_back(g_.rawSampleCT);
        g_.rawCountSNPs.push_back(p.rawCountSNP);
    }

    // Compute PGEN-specific buffer strides.
    const uint32_t keepCT = g_.keepSampleCT;
    g_.pgenGenoPtrSize =
        static_cast<int>((PgenReader::GetGenoBufPtrSize(keepCT) + 63) / 64 * 64);
    g_.pgenDosageMainPtrSize =
        static_cast<int>((PgenReader::GetDosageMainSize(keepCT) + 63) / 64 * 64);
    g_.pgenDosagePresentPtrSize =
        static_cast<int>((PgenReader::GetDosagePresentSize(keepCT) + 63) / 64 * 64);

    g_.pgenGenoBuf1PtrSize =
        static_cast<int>(
            (static_cast<std::size_t>(g_.pgenGenoPtrSize)
             + g_.pgenDosageMainPtrSize
             + g_.pgenDosagePresentPtrSize + 1 + 63) / 64 * 64);
    // No asyncBuf64 created — the coroutine frame is the buffer.
}

// --------------------------------------------------------------------------
// NOTE: callback is taken by value so it lives in the coroutine frame and
// remains valid across every co_await suspension point.  Do NOT change to a
// reference — it would dangle the moment the caller's scope is suspended.
exec::task<void> PgenBackend::stream(GenoScheduler               io_sched,
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
    const int      stride      = g_.pgenGenoBuf1PtrSize;
    const int      markerBlock = g_.numMarkerBlock;

    uint32_t   finishedMarker = 0;
    int        preFileIndex   = -1;
    int        fileIndex      = 0;
    bool       chr_ends       = false;
    uint8_t    isSexXY        = 0;
    PgenReader reader;

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
            int al_idx    = g_.marker->isEffecRevRaw(rawIndex) ? 0 : 1;
            reader.ReadDosage(
                block.buf.data() + static_cast<std::size_t>(i) * stride,
                lag_index, al_idx);
        }

        co_await stdexec::schedule(cpu_sched);
        callback(block);
        co_await stdexec::schedule(io_sched);

        finishedMarker += nextSize;
    }
} catch (...) {
    // Re-propagate through exec::task's error channel so that
    // stdexec::sync_wait re-throws on the loopDouble call-site.
    throw;
}

// --------------------------------------------------------------------------
std::unique_ptr<GenoBackend> makePgenBackend(Geno &geno)
{
    return std::make_unique<PgenBackend>(geno);
}
