#ifndef GCTA_CPU_H
#define GCTA_CPU_H

#include <vector>

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

// Portable dpotrf: Cholesky factorization of a symmetric positive-definite matrix
// (lower triangle, column-major).  Returns 0 on success, non-zero on failure.
inline int gcta_dpotrf(gcta_blas_int n, double* a, gcta_blas_int lda) {
#if defined(GCTA_USE_ACCELERATE)
    char uplo = 'L';
    gcta_blas_int info = 0;
    dpotrf_(&uplo, &n, a, &lda, &info);
    return static_cast<int>(info);
#else
    return static_cast<int>(LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', n, a, lda));
#endif
}

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

// Portable dgetrf: LU factorization with partial pivoting.
// ipiv must be pre-allocated to length n.  Returns 0 on success.
inline int gcta_dgetrf(gcta_blas_int n, double* a, gcta_blas_int lda, gcta_blas_int* ipiv) {
#if defined(GCTA_USE_ACCELERATE)
    gcta_blas_int info = 0;
    dgetrf_(&n, &n, a, &lda, ipiv, &info);
    return static_cast<int>(info);
#else
    return static_cast<int>(LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, a, lda, ipiv));
#endif
}

// Portable dgetri: in-place inversion from an LU factorization produced by gcta_dgetrf.
// ipiv must be the pivot array returned by gcta_dgetrf.  Returns 0 on success.
inline int gcta_dgetri(gcta_blas_int n, double* a, gcta_blas_int lda, gcta_blas_int* ipiv) {
#if defined(GCTA_USE_ACCELERATE)
    gcta_blas_int info = 0;
    // Workspace query.
    gcta_blas_int lwork = -1;
    double work_query = 0.0;
    dgetri_(&n, a, &lda, ipiv, &work_query, &lwork, &info);
    lwork = static_cast<gcta_blas_int>(work_query);
    std::vector<double> work(lwork);
    dgetri_(&n, a, &lda, ipiv, work.data(), &lwork, &info);
    return static_cast<int>(info);
#else
    return static_cast<int>(LAPACKE_dgetri(LAPACK_COL_MAJOR, n, a, lda, ipiv));
#endif
}

#endif  //END GCTA_CPU_H
