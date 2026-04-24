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

#elif defined(GCTA_USE_ACCELERATE)
  #ifndef EIGEN_USE_BLAS
  #define EIGEN_USE_BLAS
  #endif
  #include <cblas_new.h>
  #include <lapack.h>
  typedef __LAPACK_int gcta_blas_int;

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

  #endif
#endif

#endif  //END GCTA_CPU_H
