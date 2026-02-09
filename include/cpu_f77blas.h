#ifndef GCTA_CPU_F77BLAS_H
#define GCTA_CPU_F77BLAS_H

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
  extern "C" {
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
                const int *k, const double *alpha, const double *a, const int *lda,
                const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
    void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
                const double *a, const int *lda, const double *x, const int *incx,
                const double *beta, double *y, const int *incy);
    void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k,
                const double *alpha, const double *a, const int *lda,
                const double *beta, double *c, const int *ldc);
  }
#else
  #include <f77blas.h>
#endif

#endif  //END GCTA_F77BLAS_CPU_H