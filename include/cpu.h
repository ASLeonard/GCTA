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

#if defined(GCTA_USE_MKL)
  #ifndef EIGEN_USE_MKL_ALL
  #define EIGEN_USE_MKL_ALL
  #endif
  #include <mkl.h>
  typedef int gcta_blas_int;
#elif defined(GCTA_USE_OPENBLAS)
  #ifndef EIGEN_USE_BLAS
  #define EIGEN_USE_BLAS
  #endif
  #include <cblas.h>
  #include <lapacke.h>
  typedef lapack_int gcta_blas_int;

  // Functions with char* params need the FORTRAN_STRLEN wrapper from lapack.h
  #define dpotrf(...) LAPACK_dpotrf(__VA_ARGS__)
  #define dpotri(...) LAPACK_dpotri(__VA_ARGS__)
  #define dormqr(...) LAPACK_dormqr(__VA_ARGS__)
  // Functions without char* params map directly to the Fortran symbols
  #define dgetrf LAPACK_dgetrf
  #define dgetri LAPACK_dgetri
  #define dgeqrf LAPACK_dgeqrf
#elif defined(GCTA_USE_ACCELERATE)
  #ifndef EIGEN_USE_BLAS
  #define EIGEN_USE_BLAS
  #endif
  #include <cblas_new.h>
  #include <lapack.h>
  typedef __LAPACK_int gcta_blas_int;

  #define dpotrf dpotrf_
  #define dpotri dpotri_
  #define dgetrf dgetrf_
  #define dgetri dgetri_
  #define dgeqrf dgeqrf_
  #define dormqr dormqr_
#else
  #ifndef EIGEN_USE_BLAS
  #define EIGEN_USE_BLAS
  #endif

  #if defined(GCTA_USE_ACCELERATE)
    #include <cblas_new.h>
    #include <lapack.h>
    typedef __LAPACK_int gcta_blas_int;
  #else
    #include <cblas.h>
    #if defined(__has_include)
      #if __has_include(<lapacke.h>)
        #include <lapacke.h>
      #endif
    #endif
    typedef int gcta_blas_int;

    // Fortran LAPACK declarations for generic BLAS/LAPACK
    extern "C" {
      void dpotrf_(const char *uplo, const gcta_blas_int *n, double *a,
                   const gcta_blas_int *lda, gcta_blas_int *info);
      void dpotri_(const char *uplo, const gcta_blas_int *n, double *a,
                   const gcta_blas_int *lda, gcta_blas_int *info);
      void dgetrf_(const gcta_blas_int *m, const gcta_blas_int *n, double *a,
                   const gcta_blas_int *lda, gcta_blas_int *ipiv, gcta_blas_int *info);
      void dgetri_(const gcta_blas_int *n, double *a, const gcta_blas_int *lda,
                   const gcta_blas_int *ipiv, double *work, const gcta_blas_int *lwork, gcta_blas_int *info);
      void dgeqrf_(const gcta_blas_int *m, const gcta_blas_int *n, double *a,
                   const gcta_blas_int *lda, double *tau, double *work, const gcta_blas_int *lwork, gcta_blas_int *info);
      void dormqr_(const char *side, const char *trans, const gcta_blas_int *m,
                   const gcta_blas_int *n, const gcta_blas_int *k, const double *a,
                   const gcta_blas_int *lda, const double *tau, double *c, const gcta_blas_int *ldc,
                   double *work, const gcta_blas_int *lwork, gcta_blas_int *info);
    }
  #endif

  // Alias underscore-suffixed Fortran names to the bare names used in this codebase.
  // Must come after includes so Accelerate's <lapack.h> can set up its own ILP64 aliases first.
  #define dpotrf dpotrf_
  #define dpotri dpotri_
  #define dgetrf dgetrf_
  #define dgetri dgetri_
  #define dgeqrf dgeqrf_
  #define dormqr dormqr_
#endif

#endif  //END GCTA_CPU_H
