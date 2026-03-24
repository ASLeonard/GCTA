#ifndef GCTA_CPU_F77BLAS_H
#define GCTA_CPU_F77BLAS_H

#include "cpu.h"

#if defined(__x86_64__) || (defined(_M_X64) && !defined(_M_ARM64EC)) || defined(__amd64)
  #define GCTA_ARCH_x86_64 1
#else
  #define GCTA_ARCH_x86_64 0
#endif

#if defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__i386)
  #define GCTA_ARCH_i386 1
#else
  #define GCTA_ARCH_i386 0
#endif

#if GCTA_ARCH_x86_64 || GCTA_ARCH_i386
  #define GCTA_CPU_x86 1
#else
  #define GCTA_CPU_x86 0
#endif

#if defined(__arm__)
  #define GCTA_ARCH_ARM 1
#else
  #define GCTA_ARCH_ARM 0
#endif

#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
  #define GCTA_ARCH_ARM64 1
#else
  #define GCTA_ARCH_ARM64 0
#endif

#if GCTA_ARCH_ARM || GCTA_ARCH_ARM64
  #define GCTA_CPU_ARM 1
#else
  #define GCTA_CPU_ARM 0
#endif

#if defined(GCTA_USE_MKL)
  #ifndef EIGEN_USE_MKL_ALL
  #define EIGEN_USE_MKL_ALL
  #endif
  #include <mkl.h>
#elif defined(GCTA_USE_OPENBLAS)
  #include <f77blas.h>
#elif defined(GCTA_USE_ACCELERATE)
  static inline enum CBLAS_TRANSPOSE gcta_cblas_transpose(char trans) {
    return (trans == 'T' || trans == 't') ? CblasTrans : CblasNoTrans;
  }

  static inline enum CBLAS_UPLO gcta_cblas_uplo(char uplo) {
    return (uplo == 'U' || uplo == 'u') ? CblasUpper : CblasLower;
  }

  static inline void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
                            const int *k, const double *alpha, const double *a, const int *lda,
                            const double *b, const int *ldb, const double *beta, double *c, const int *ldc) {
    const gcta_blas_int m64 = (gcta_blas_int)*m;
    const gcta_blas_int n64 = (gcta_blas_int)*n;
    const gcta_blas_int k64 = (gcta_blas_int)*k;
    const gcta_blas_int lda64 = (gcta_blas_int)*lda;
    const gcta_blas_int ldb64 = (gcta_blas_int)*ldb;
    const gcta_blas_int ldc64 = (gcta_blas_int)*ldc;
    cblas_dgemm(CblasColMajor, gcta_cblas_transpose(*transa), gcta_cblas_transpose(*transb),
                m64, n64, k64, *alpha, a, lda64, b, ldb64, *beta, c, ldc64);
  }

  static inline void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
                            const double *a, const int *lda, const double *x, const int *incx,
                            const double *beta, double *y, const int *incy) {
    const gcta_blas_int m64 = (gcta_blas_int)*m;
    const gcta_blas_int n64 = (gcta_blas_int)*n;
    const gcta_blas_int lda64 = (gcta_blas_int)*lda;
    const gcta_blas_int incx64 = (gcta_blas_int)*incx;
    const gcta_blas_int incy64 = (gcta_blas_int)*incy;
    cblas_dgemv(CblasColMajor, gcta_cblas_transpose(*trans), m64, n64, *alpha, a, lda64, x, incx64, *beta, y, incy64);
  }

  static inline void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k,
                            const double *alpha, const double *a, const int *lda,
                            const double *beta, double *c, const int *ldc) {
    const gcta_blas_int n64 = (gcta_blas_int)*n;
    const gcta_blas_int k64 = (gcta_blas_int)*k;
    const gcta_blas_int lda64 = (gcta_blas_int)*lda;
    const gcta_blas_int ldc64 = (gcta_blas_int)*ldc;
    cblas_dsyrk(CblasColMajor, gcta_cblas_uplo(*uplo), gcta_cblas_transpose(*trans),
                n64, k64, *alpha, a, lda64, *beta, c, ldc64);
  }
#else
  #include <f77blas.h>
#endif

#endif  //END GCTA_F77BLAS_CPU_H