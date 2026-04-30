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

#elif defined(GCTA_USE_ACCELERATE)
  #include <cblas_new.h>
  #include <lapack.h>
  typedef __LAPACK_int gcta_blas_int;
#endif

#endif  //END GCTA_CPU_H
