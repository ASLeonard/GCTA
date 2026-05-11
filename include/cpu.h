#ifndef GCTA_CPU_H
#define GCTA_CPU_H

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64EC)) || defined(__amd64)
  #define GCTA_CPU_x86 1
#else
  #define GCTA_CPU_x86 0
#endif

#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
  #define GCTA_CPU_ARM 1
#else
  #define GCTA_CPU_ARM 0
#endif

//Can replace this header entirely with a CMake file(GENERATE) approach:

#if defined(GCTA_USE_MKL)
  #include <mkl.h>
  typedef int gcta_blas_int;
#elif defined(GCTA_USE_OPENBLAS)
  #include <cblas.h>
  #include <lapacke.h>
  typedef lapack_int gcta_blas_int;

#elif defined(GCTA_USE_AOCL)
  #include <cblas.h>
  #include <lapacke.h>
  typedef lapack_int gcta_blas_int;

#elif defined(GCTA_USE_ACCELERATE)
  #include <cblas_new.h>
  #include <lapack.h>
  typedef __LAPACK_int gcta_blas_int;
#endif

// Portable dpotri: in-place inversion of a Cholesky-factored symmetric positive-definite
// matrix (lower triangle).  Returns 0 on success, non-zero on failure.
inline int gcta_dpotri(gcta_blas_int n, double* a, gcta_blas_int lda) {
#if defined(GCTA_USE_ACCELERATE)
    // Accelerate exposes Fortran-ABI dpotri_ with pointer arguments.
    char uplo = 'L';
    gcta_blas_int info = 0;
    dpotri_(&uplo, &n, a, &lda, &info);
    return static_cast<int>(info);
#else
    // MKL, OpenBLAS, and AOCL all provide the LAPACKE C interface.
    return static_cast<int>(LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', n, a, lda));
#endif
}

#endif  //END GCTA_CPU_H
