/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Mocks the same memory allocation api cross platform

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GCTA2_MEM_HPP
#define GCTA2_MEM_HPP

#ifdef _WIN32
#include <malloc.h>
#define posix_memalign(p, a, s) ( ((*(p)) = _aligned_malloc((s), (a))), *(p) ? 0 : errno )
#define posix_mem_free _aligned_free
#else
#include <cstdio>
#include <stdlib.h>
#define posix_mem_free free
#endif

#include <cstddef>
#include <new>

// ---------------------------------------------------------------------------
// AlignedAllocator<T, Align>
//
// Drop-in std::allocator replacement that guarantees Align-byte alignment.
// Intended for std::vector buffers passed to plink2/SIMD code that emits
// vmovdqa (aligned move) on Linux/GCC and would segfault on misaligned data.
// ---------------------------------------------------------------------------
template<typename T, std::size_t Align = 32>
struct AlignedAllocator {
    using value_type = T;

    AlignedAllocator() noexcept = default;
    template<typename U>
    AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}

    T* allocate(std::size_t n) {
        if (n == 0) return nullptr;
        void* p = nullptr;
        if (posix_memalign(&p, Align, n * sizeof(T)) != 0)
            throw std::bad_alloc();
        return static_cast<T*>(p);
    }

    void deallocate(T* p, std::size_t) noexcept {
        posix_mem_free(p);
    }

    template<typename U>
    struct rebind { using other = AlignedAllocator<U, Align>; };
};

template<typename T, typename U, std::size_t A>
bool operator==(const AlignedAllocator<T,A>&, const AlignedAllocator<U,A>&) noexcept { return true; }
template<typename T, typename U, std::size_t A>
bool operator!=(const AlignedAllocator<T,A>&, const AlignedAllocator<U,A>&) noexcept { return false; }

// These functions are only for test purpose, don't forget to remove calls
int getVMemKB();
int getMemKB();

int getVMPeakKB();
int getMemPeakKB();

#endif //GCTA2_MEM_HPP


