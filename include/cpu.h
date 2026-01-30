#ifndef GCTA_CPU_H
#define GCTA_CPU_H

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

// Try to detect if we're using MKL or generic BLAS/LAPACK
// MKL is preferred on x86 if available, but not required
#if GCTA_CPU_x86
  // Check if MKL is available (set by CMake via EIGEN_USE_MKL_ALL)
  #ifdef EIGEN_USE_MKL_ALL
    #include <mkl.h>
  #else
    // Fall back to generic BLAS/LAPACK even on x86
    #ifndef EIGEN_USE_BLAS
      #define EIGEN_USE_BLAS
    #endif
    #include <cblas.h>
    #if defined(__APPLE__)
      #include <clapack.h>
      // On macOS Accelerate, BLAS Fortran functions don't have trailing underscore
      extern "C" {
        void dsyrk(const char* uplo, const char* trans, const int* n, const int* k,
                   const double* alpha, const double* a, const int* lda,
                   const double* beta, double* c, const int* ldc);
        void dgemm(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                   const double* alpha, const double* a, const int* lda,
                   const double* b, const int* ldb,
                   const double* beta, double* c, const int* ldc);
      }
      // Don't define macros - use the functions directly
    #else
      extern "C" {
        // BLAS function declarations (with trailing underscore on Linux)
        void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
                   const double* alpha, const double* a, const int* lda,
                   const double* beta, double* c, const int* ldc);
        void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                   const double* alpha, const double* a, const int* lda,
                   const double* b, const int* ldb,
                   const double* beta, double* c, const int* ldc);
        
        // LAPACK function declarations
        void dgeqrf_(const int* m, const int* n, double* a, const int* lda,
                    double* tau, double* work, const int* lwork, int* info);
        void dormqr_(const char* side, const char* trans, const int* m, const int* n,
                    const int* k, const double* a, const int* lda, const double* tau,
                    double* c, const int* ldc, double* work, const int* lwork, int* info);
      }
      // Define macros to map non-underscore names to underscore names
      #define dsyrk dsyrk_
      #define dgemm dgemm_
      #define dgeqrf dgeqrf_
      #define dormqr dormqr_
    #endif
  #endif
#else
  // ARM and other architectures use generic BLAS/LAPACK
  #ifndef EIGEN_USE_BLAS
    #define EIGEN_USE_BLAS
  #endif
  #include <cblas.h>
  #if defined(__APPLE__)
    // macOS provides clapack.h with Accelerate framework
    #include <clapack.h>
    // On macOS Accelerate, BLAS Fortran functions don't have trailing underscore
    extern "C" {
      void dsyrk(const char* uplo, const char* trans, const int* n, const int* k,
                 const double* alpha, const double* a, const int* lda,
                 const double* beta, double* c, const int* ldc);
      void dgemm(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                 const double* alpha, const double* a, const int* lda,
                 const double* b, const int* ldb,
                 const double* beta, double* c, const int* ldc);
    }
    // Don't define macros - use the functions directly
  #else
    // Generic BLAS/LAPACK
    extern "C" {
      // BLAS function declarations (with trailing underscore on Linux)
      void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
                 const double* alpha, const double* a, const int* lda,
                 const double* beta, double* c, const int* ldc);
      void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
                 const double* alpha, const double* a, const int* lda,
                 const double* b, const int* ldb,
                 const double* beta, double* c, const int* ldc);
      
      // LAPACK function declarations
      void dgeqrf_(const int* m, const int* n, double* a, const int* lda,
                  double* tau, double* work, const int* lwork, int* info);
      void dormqr_(const char* side, const char* trans, const int* m, const int* n,
                  const int* k, const double* a, const int* lda, const double* tau,
                  double* c, const int* ldc, double* work, const int* lwork, int* info);
    }
    // Define macros to map non-underscore names to underscore names
    #define dsyrk dsyrk_
    #define dgemm dgemm_
    #define dgeqrf dgeqrf_
    #define dormqr dormqr_
  #endif
#endif

#endif  //END GCTA_CPU_H
