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
#include <array>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <mutex>
#include <numeric>
#include <span>
#include <thread>
#include <vector>

// --------------------------------------------------------------------------
// Helper: aligned-allocated mask buffer with posix_mem_free deleter.
// --------------------------------------------------------------------------
namespace bed_detail {
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
} // namespace bed_detail

using bed_detail::MaskBuf;
using bed_detail::allocMask;

// --------------------------------------------------------------------------
// BedBackend
// --------------------------------------------------------------------------
class BedBackend : public GenoBackend {
public:
    explicit BedBackend(Geno &geno);
    ~BedBackend() override;

    exec::task<void> stream(GenoScheduler               io_sched,
                            GenoScheduler               cpu_sched,
                            std::vector<uint32_t>       extractIndex,
                            Callback                    callback) override;

private:
    Geno &g_;

    // Owned mask buffers (lifetime = this backend).
    // The raw pointers in Geno (g_.keepMaskPtr etc.) point into these.
    MaskBuf keepMask_;
    MaskBuf keepMaskInter_;
    MaskBuf sexMask_;
    MaskBuf sexMaskInter_;
    MaskBuf heterogameticMask_;
    MaskBuf heterogameticMaskInter_;

    // Two reusable slots keep allocator pressure off the hot loop.
    std::array<GenoBlock, 2> blocks_;

    // stream() is not safe to run concurrently on the same backend instance.
    std::atomic<bool> streamActive_{false};
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

    const std::size_t maxBlockMarkers =
        static_cast<std::size_t>(std::max(1, g_.numMarkerBlock));
    const std::size_t blockElemCapacity =
        static_cast<std::size_t>(g_.bedRawGenoBuf1PtrSize) * maxBlockMarkers;
    for (auto &block : blocks_) {
        block.buf.resize(blockElemCapacity);
    }

    std::size_t nMask = static_cast<std::size_t>(g_.maskPtrSize);
    keepMask_       = allocMask(nMask);
    keepMaskInter_  = allocMask(nMask);
    sexMask_        = allocMask(nMask);
    sexMaskInter_   = allocMask(nMask);
    heterogameticMask_       = allocMask(nMask);
    heterogameticMaskInter_  = allocMask(nMask);

    if (!keepMask_ || !keepMaskInter_ ||
        !sexMask_  || !sexMaskInter_  ||
        !heterogameticMask_ || !heterogameticMaskInter_) {
        LOGGER.e(0, "BedBackend: can't allocate mask buffers.");
    }

    // Expose raw pointers so that Geno's getGenoDouble_bed can use them.
    g_.keepMaskPtr      = keepMask_.get();
    g_.keepMaskInterPtr = keepMaskInter_.get();
    g_.sexMaskPtr       = sexMask_.get();
    g_.sexMaskInterPtr  = sexMaskInter_.get();
    g_.heterogameticMaskPtr      = heterogameticMask_.get();
    g_.heterogameticMaskInterPtr = heterogameticMaskInter_.get();

    PgenReader::SetSampleSubsets(g_.sampleKeepIndex, rawCT,
                                 g_.keepMaskPtr, g_.keepMaskInterPtr);
    PgenReader::SetSampleSubsets(g_.keepSexIndex,  rawCT,
                                 g_.sexMaskPtr,  g_.sexMaskInterPtr);
    PgenReader::SetSampleSubsets(g_.keepHeterogameticIndex, rawCT,
                                 g_.heterogameticMaskPtr, g_.heterogameticMaskInterPtr);
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
    g_.heterogameticMaskPtr      = nullptr;
    g_.heterogameticMaskInterPtr = nullptr;
}

// --------------------------------------------------------------------------
// NOTE: callback is taken by value so it lives in the coroutine frame and
// remains valid across every co_await suspension point.  Do NOT change to a
// reference — it would dangle the moment the caller's scope is suspended.
exec::task<void> BedBackend::stream(GenoScheduler               io_sched,
                                    GenoScheduler               cpu_sched,
                                    std::vector<uint32_t>       extractIndex,
                                    Callback                    callback)
try {
    (void)cpu_sched;

    if (streamActive_.exchange(true, std::memory_order_acq_rel)) {
        LOGGER.e(0, "BedBackend::stream called concurrently on one backend instance.");
    }
    struct StreamActiveReset {
        std::atomic<bool> &active;
        ~StreamActiveReset() { active.store(false, std::memory_order_release); }
    } streamActiveReset{streamActive_};

    // ── Transfer to the I/O thread ───────────────────────────────────────
    co_await stdexec::schedule(io_sched);

    // ── Raw marker map and small per-block raw-index scratch ─────────────
    const std::vector<uint32_t> &raw_marker_index = g_.marker->get_extract_index();

    const uint32_t numMarker    = static_cast<uint32_t>(extractIndex.size());
    const int      stride       = g_.bedRawGenoBuf1PtrSize;
    const int      markerBlock  = g_.numMarkerBlock;
    std::vector<uint32_t> rawRefScratch;
    rawRefScratch.reserve(static_cast<std::size_t>(std::max(1, markerBlock)));
    std::vector<int> lagIndexScratch;
    lagIndexScratch.reserve(static_cast<std::size_t>(std::max(1, markerBlock)));

    uint32_t finishedMarker = 0;
    int      preFileIndex   = -1;
    int      fileIndex      = 0;
    bool     chr_ends       = false;
    uint8_t  sexChromType        = 0;
    PgenReader reader;

    std::mutex queueMutex;
    std::condition_variable_any readyCv;
    std::condition_variable_any freeCv;
    std::deque<int> readySlots;
    std::array<bool, 2> slotInUse{false, false};
    bool done = false;
    std::exception_ptr callbackError;
    std::stop_source stopSource;
    const std::stop_token stopToken = stopSource.get_token();

    std::jthread computeThread([&]() {
        while (true) {
            int slot = -1;
            {
                std::unique_lock<std::mutex> lock(queueMutex);
                readyCv.wait(lock, stopToken, [&]() {
                    return callbackError || done || !readySlots.empty();
                });

                if (stopToken.stop_requested()) break;
                if (callbackError) break;
                if (readySlots.empty()) {
                    if (done) break;
                    continue;
                }

                slot = readySlots.front();
                readySlots.pop_front();
            }

            try {
                callback(blocks_[static_cast<std::size_t>(slot)]);
            } catch (...) {
                std::lock_guard<std::mutex> lock(queueMutex);
                callbackError = std::current_exception();
                done = true;
                stopSource.request_stop();
                readyCv.notify_all();
                freeCv.notify_all();
                break;
            }

            {
                std::lock_guard<std::mutex> lock(queueMutex);
                slotInUse[static_cast<std::size_t>(slot)] = false;
            }
            freeCv.notify_one();
        }
    });

    // ── Main I/O loop ────────────────────────────────────────────────────
    while (finishedMarker < numMarker) {
        rawRefScratch.clear();
        const uint32_t remain = numMarker - finishedMarker;
        const uint32_t probeSize = std::min<uint32_t>(
            static_cast<uint32_t>(std::max(1, markerBlock)), remain);
        rawRefScratch.resize(static_cast<std::size_t>(probeSize));
        for (uint32_t i = 0; i < probeSize; ++i) {
            rawRefScratch[static_cast<std::size_t>(i)] =
                raw_marker_index[extractIndex[finishedMarker + i]];
        }

        uint32_t nextSize = g_.marker->getNextSize(
            rawRefScratch, 0, probeSize,
            fileIndex, chr_ends, sexChromType);
        if (nextSize == 0) break;

        int slot = -1;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            freeCv.wait(lock, stopToken, [&]() {
                return callbackError || !slotInUse[0] || !slotInUse[1];
            });
            if (stopToken.stop_requested() && callbackError) {
                std::rethrow_exception(callbackError);
            }
            if (callbackError) std::rethrow_exception(callbackError);

            slot = !slotInUse[0] ? 0 : 1;
            slotInUse[static_cast<std::size_t>(slot)] = true;
        }

        // Fill a reusable block on the I/O thread.
        GenoBlock &block = blocks_[static_cast<std::size_t>(slot)];
        block.numMarkers = nextSize;
        block.sexChromType = sexChromType;
        block.fileIndex = fileIndex;
        const std::size_t requiredElems =
            static_cast<std::size_t>(stride) * nextSize;
        if (block.buf.size() < requiredElems) {
            block.buf.resize(requiredElems);
        }

        block.extractIndex = std::span<const uint32_t>{
            extractIndex.data() + finishedMarker,
            static_cast<std::size_t>(nextSize)};

        if (preFileIndex != fileIndex) {
            reader.Load(g_.geno_files[fileIndex],
                        &g_.rawCountSamples[fileIndex],
                        &g_.rawCountSNPs[fileIndex],
                        g_.sampleKeepIndex);
            preFileIndex = fileIndex;
        }

        lagIndexScratch.resize(static_cast<std::size_t>(nextSize));
        const int baseIndex = g_.baseIndexLookup[fileIndex];
        bool contiguousLag = (nextSize > 0);
        for (uint32_t i = 0; i < nextSize; ++i) {
            lagIndexScratch[static_cast<std::size_t>(i)] =
                static_cast<int>(rawRefScratch[static_cast<std::size_t>(i)]) - baseIndex;
            if (i > 0 && contiguousLag) {
                contiguousLag =
                    (lagIndexScratch[static_cast<std::size_t>(i)]
                     == lagIndexScratch[static_cast<std::size_t>(i - 1)] + 1);
            }
        }
        if (contiguousLag) {
            reader.ReadRawFullHardRange(
                block.buf.data(), lagIndexScratch.front(), static_cast<int>(nextSize), stride);
        } else {
            reader.ReadRawFullHardBatch(
                block.buf.data(), lagIndexScratch.data(), static_cast<int>(nextSize), stride);
        }

        {
            std::lock_guard<std::mutex> lock(queueMutex);
            readySlots.push_back(slot);
        }
        readyCv.notify_one();

        finishedMarker += nextSize;
    }

    {
        std::lock_guard<std::mutex> lock(queueMutex);
        done = true;
    }
    readyCv.notify_all();

    if (computeThread.joinable()) {
        computeThread.join();
    }
    if (callbackError) {
        std::rethrow_exception(callbackError);
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
