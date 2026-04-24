/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Coroutine-based genotype I/O backend interface.

   Replaces the six dispatch-map families (preGenoDoubleFuncs,
   getGenoDoubleFuncs, endGenoDoubleFuncs, readGenoFuncs, preProcessFuncs,
   endProcessFuncs) and the AsyncBuffer<uintptr_t> pipeline with a
   single virtual stream() method per format.

   Each backend:
     - runs pre-logic (mask allocation etc.) in its constructor
     - runs end-logic (mask deallocation etc.) in its destructor
     - implements stream() as a C++ coroutine that reads raw genotype blocks
       on an I/O thread and hands them to the user callback on a CPU thread,
       using the NVIDIA/stdexec P2300 reference implementation.

   The callback retains the classic  void(uintptr_t* buf, const std::vector<uint32_t>& exIndex)
   contract — loopDouble wraps it so that no callers outside Geno need to change.
*/
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

// --------------------------------------------------------------------------
// stdexec (P2300 reference implementation) — required for exec::task<void>
// and exec::static_thread_pool.
// --------------------------------------------------------------------------
#include <stdexec/execution.hpp>
#include <exec/task.hpp>
#include <exec/static_thread_pool.hpp>

// --------------------------------------------------------------------------
// GenoBlock
//
// One block of raw genotype data as filled by a backend's I/O coroutine.
// The layout of buf[] is format-specific (BED / BGEN / PGEN); callers that
// need to interpret it call Geno::getGenoDouble(buf.data(), i, &item).
// --------------------------------------------------------------------------
struct GenoBlock {
    std::vector<uintptr_t> buf;          ///< raw packed data, stride × numMarkers elements
    std::vector<uint32_t>  extractIndex; ///< marker extract-indices for this block
    uint32_t               numMarkers = 0;
    uint8_t                isSexXY   = 0; ///< 0=autosome, 1=chrX, 2=chrY
    int                    fileIndex  = 0; ///< which geno_file[] this block came from
};

// --------------------------------------------------------------------------
// GenoScheduler — the concrete scheduler type from exec::static_thread_pool.
// Using decltype avoids hard-coding an implementation-internal type name.
// --------------------------------------------------------------------------
using GenoScheduler =
    decltype(std::declval<exec::static_thread_pool &>().get_scheduler());

// --------------------------------------------------------------------------
// GenoBackend — abstract base for BED / BGEN / PGEN backends.
//
// Lifecycle:
//   construction  → runs format-specific pre-logic  (mask allocation, file setup)
//   stream()      → coroutine I/O loop
//   destruction   → runs format-specific end-logic  (mask deallocation, file close)
//
// All three backends share this interface so that loopDouble() can hold a
// std::unique_ptr<GenoBackend> without knowing the concrete format.
// --------------------------------------------------------------------------
class GenoBackend {
public:
    using Callback = std::function<void(GenoBlock &)>;

    GenoBackend()          = default;
    virtual ~GenoBackend() = default;

    // Non-copyable, non-moveable (owns OS resources and Geno pointers).
    GenoBackend(const GenoBackend &)            = delete;
    GenoBackend &operator=(const GenoBackend &) = delete;

    /// Stream all blocks for extractIndex through callback.
    /// I/O runs on io_sched (1 thread); callback runs on cpu_sched.
    /// The coroutine returns only when all blocks have been processed.
    virtual exec::task<void> stream(
        GenoScheduler               io_sched,
        GenoScheduler               cpu_sched,
        const std::vector<uint32_t> &extractIndex,
        Callback                    callback) = 0;
};
