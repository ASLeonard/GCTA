CMAKE_MINIMUM_REQUIRED(VERSION 3.28)

project(GCTA2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 23)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_COMPILER_IS_GNUCC)
    message("Compiling by GCC ${CMAKE_CXX_COMPILER_VERSION}")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
        message(FATAL_ERROR "GCTA requires GCC >= 13 for reliable C++23 support (found ${CMAKE_CXX_COMPILER_VERSION})")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiling by Clang ${CMAKE_CXX_COMPILER_VERSION}")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0)
        message(FATAL_ERROR "GCTA requires Clang >= 17 for reliable C++23 support (found ${CMAKE_CXX_COMPILER_VERSION})")
    endif()
    if(APPLE AND CMAKE_CXX_COMPILER MATCHES "/opt/homebrew/opt/llvm")
        set(_llvm_prefix "/opt/homebrew/opt/llvm")
        string(APPEND CMAKE_CXX_FLAGS " -stdlib=libc++")
        string(APPEND CMAKE_EXE_LINKER_FLAGS " -L${_llvm_prefix}/lib/c++")
        link_directories("${_llvm_prefix}/lib/c++")
        list(APPEND CMAKE_BUILD_RPATH "${_llvm_prefix}/lib/c++")
        message(STATUS "Using Homebrew LLVM libc++ from ${_llvm_prefix}/lib/c++")
    endif()
else()
    message(FATAL_ERROR "GCTA can only compile with GCC or Clang")
endif()

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
    message(STATUS "Defaulting to Release build")
endif()

# Native arch optimisations: Apple Silicon uses -mcpu, x86_64 uses -march
if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease$")
    set(_arch_flag $<IF:$<BOOL:${APPLE}>,-mcpu=native,-march=native>)
    add_compile_options(-O3 ${_arch_flag} -pipe)
endif()

set(GCTA_BLAS_BACKEND "Auto" CACHE STRING "BLAS/LAPACK backend: Auto, MKL, OpenBLAS, or Accelerate")
set_property(CACHE GCTA_BLAS_BACKEND PROPERTY STRINGS Auto MKL OpenBLAS Accelerate)
option(GCTA_LIST_TARGETS "List imported and local targets during configure" OFF)
option(GCTA_UNITY_BUILD "Enable unity builds for faster compilation" ON)
option(GCTA_ENABLE_IPO "Enable IPO/LTO in Release builds" ON)

# Initialize plink-ng submodule if needed
if(EXISTS "${CMAKE_SOURCE_DIR}/.git" AND NOT EXISTS "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include/pgenlib_ffi_support.h")
    find_package(Git QUIET)
    if(GIT_FOUND)
        message(STATUS "Initializing plink-ng submodule...")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}")
        endif()
    else()
        message(FATAL_ERROR "plink-ng submodule not initialized. Please run: git submodule update --init --recursive")
    endif()
endif()

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include/pgenlib_ffi_support.h")
    message(FATAL_ERROR "plink-ng submodule is missing. Please run: git submodule update --init --recursive")
endif()

# Fetch 3rd parties
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
    zlib
    URL https://github.com/madler/zlib/archive/refs/tags/v1.3.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL
    FIND_PACKAGE_ARGS NAMES ZLIB)
FetchContent_Declare(
    eigen
    URL https://github.com/eigen-mirror/eigen/archive/refs/tags/5.0.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(
    spectra
    URL https://github.com/yixuan/spectra/archive/refs/tags/v1.2.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(
    zstd
    URL https://github.com/facebook/zstd/archive/refs/tags/v1.5.7.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR "build/cmake"
    EXCLUDE_FROM_ALL
    FIND_PACKAGE_ARGS NAMES zstd)
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-cmake.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL)
FetchContent_Declare(
    sqlite
    URL https://www.sqlite.org/2026/sqlite-amalgamation-3510200.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL)

# zlib must be fetched before Boost so Boost iostreams can find it at configure time.
# Kept outside block() so zlib_POPULATED and zlib_SOURCE_DIR are visible in parent scope.
set(ZLIB_BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(zlib)

# After zlib is available, set hints so Boost iostreams links against it correctly.
if(zlib_POPULATED)
    set(ZLIB_INCLUDE_DIR "${zlib_SOURCE_DIR}" CACHE PATH "" FORCE)
    if(UNIX)
        set(ZLIB_LIBRARY "${zlib_BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}z${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE FILEPATH "" FORCE)
    else()
        set(ZLIB_LIBRARY "${zlib_BINARY_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}zlibstatic${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE FILEPATH "" FORCE)
    endif()
endif()

block()
    set(ZSTD_BUILD_PROGRAMS OFF)
    set(ZSTD_BUILD_TESTS OFF)
    set(ZSTD_BUILD_SHARED ON)
    set(BOOST_ENABLE_CMAKE ON)
    set(BOOST_INCLUDE_LIBRARIES "algorithm;math;crc;iostreams")
    set(SPECTRA_BUILD_TESTS OFF)
    set(BUILD_TESTING OFF)
    # Explicitly enable zlib support in Boost iostreams and point it at our zlib.
    # Without this Boost may silently disable the zlib filter if it doesn't find
    # zlib on its own, producing the missing zlib_error::check link error later.
    set(BOOST_IOSTREAMS_ENABLE_ZLIB ON)
    if(zlib_POPULATED)
        set(ZLIB_ROOT "${zlib_BINARY_DIR}" CACHE PATH "" FORCE)
        set(ZLIB_INCLUDE_DIR "${zlib_SOURCE_DIR};${zlib_BINARY_DIR}" CACHE PATH "" FORCE)
    endif()
    FetchContent_MakeAvailable(zstd Boost)
endblock()

# Eigen and Spectra have entangled install rules (Spectra exports reference
# Eigen targets) that break with EXCLUDE_FROM_ALL. Suppress via
# CMAKE_SKIP_INSTALL_RULES instead, which must be set in the parent scope.
set(CMAKE_SKIP_INSTALL_RULES ON)
FetchContent_MakeAvailable(eigen spectra)
set(CMAKE_SKIP_INSTALL_RULES OFF)

# sqlite has no options to isolate, and sqlite_SOURCE_DIR must be visible
# in the parent scope for the add_library(sqlite3) call below.
FetchContent_MakeAvailable(sqlite)

# Ensure a consistent ZLIB::ZLIB target exists for downstream linking.
# - System zlib: find_package (via FIND_PACKAGE_ARGS) already created ZLIB::ZLIB
# - Fetched zlib: alias the real build target so the library is built before use
if(zlib_POPULATED)
    set(GCTA_ZLIB_INCLUDE_DIRS "${zlib_SOURCE_DIR}" "${zlib_BINARY_DIR}")
    if(NOT TARGET ZLIB::ZLIB)
        # zlib's CMakeLists creates 'zlib' (shared) and 'zlibstatic' on all platforms,
        # and 'zlib' (static, same name) on UNIX. Alias whichever is available.
        if(TARGET zlibstatic)
            add_library(ZLIB::ZLIB ALIAS zlibstatic)
        elseif(TARGET zlib)
            add_library(ZLIB::ZLIB ALIAS zlib)
        endif()
    endif()
else()
    # System zlib — find_package already created ZLIB::ZLIB
    set(GCTA_ZLIB_INCLUDE_DIRS "${ZLIB_INCLUDE_DIRS}")
endif()

# Resolve zstd link target — imported name differs between system and fetched.
if(zstd_POPULATED AND NOT zstd_FOUND)
    set(GCTA_ZSTD_INCLUDE_DIRS "${zstd_SOURCE_DIR}/lib")
    set(GCTA_ZSTD_LIBS libzstd)
else()
    set(GCTA_ZSTD_INCLUDE_DIRS "")
    set(GCTA_ZSTD_LIBS zstd::libzstd_shared)
endif()

add_library(sqlite3 STATIC ${sqlite_SOURCE_DIR}/sqlite3.c)
target_include_directories(sqlite3 PUBLIC ${sqlite_SOURCE_DIR})

find_package(Threads REQUIRED)

if(NOT WIN32)
    # NOTE: Install gsl with your system package manager
    find_package(GSL REQUIRED)
endif()

# Find OpenMP
if(APPLE)
    # Use manual flags for OpenMP with Apple Clang
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_C_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()

if(WIN32)
    set(OpenMP_omp_LIBRARY "${LLVM_ROOT}/lib/libomp.lib")
    set(OpenMP_C_FLAGS "-fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_FLAGS "-fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_INCLUDE_DIR "${LLVM_ROOT}/lib/clang/20/include")
endif()

find_package(OpenMP REQUIRED)

set(OPENBLAS_INCLUDE_DIR "" CACHE PATH "Path to OpenBLAS headers (used when GCTA_BLAS_BACKEND=OpenBLAS)")
set(GCTA_BLAS_LIBS "")
set(GCTA_BLAS_INCLUDE_DIRS "")
set(GCTA_BLAS_DEFINES "")

# Resolve "Auto" to a concrete backend
if(GCTA_BLAS_BACKEND STREQUAL "Auto")
    if(APPLE)
        set(_GCTA_BLAS "Accelerate")
    else()
        set(_GCTA_BLAS "MKL")
    endif()
else()
    set(_GCTA_BLAS "${GCTA_BLAS_BACKEND}")
endif()

# Validate platform constraints
if(APPLE AND NOT CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    message(FATAL_ERROR "GCTA macOS support is arm64 only.")
endif()
if(APPLE AND _GCTA_BLAS STREQUAL "MKL")
    message(FATAL_ERROR "MKL is not supported on macOS. Use Accelerate or OpenBLAS.")
endif()
if(NOT APPLE AND _GCTA_BLAS STREQUAL "Accelerate")
    message(FATAL_ERROR "Accelerate is only available on macOS. Use MKL or OpenBLAS.")
endif()

if(_GCTA_BLAS STREQUAL "Accelerate")
    find_library(ACCELERATE_LIBRARY Accelerate REQUIRED)
    list(APPEND GCTA_BLAS_LIBS ${ACCELERATE_LIBRARY})
    list(APPEND GCTA_BLAS_DEFINES GCTA_USE_ACCELERATE ACCELERATE_NEW_LAPACK EIGEN_USE_BLAS EIGEN_USE_LAPACK)
    set(GCTA_BLAS_PROVIDER "Accelerate" CACHE INTERNAL "")
    if(NOT CMAKE_OSX_SYSROOT)
        execute_process(
            COMMAND xcrun --sdk macosx --show-sdk-path
            OUTPUT_VARIABLE MACOSX_SDK_PATH
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        set(CMAKE_OSX_SYSROOT "${MACOSX_SDK_PATH}")
    endif()
    set(VECLIB_INCLUDE_DIR "${CMAKE_OSX_SYSROOT}/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers")
    list(APPEND GCTA_BLAS_INCLUDE_DIRS ${VECLIB_INCLUDE_DIR})

elseif(_GCTA_BLAS STREQUAL "OpenBLAS")
    # Check Homebrew path first (macOS), then fall back to find_package
    set(_openblas_prefix "/opt/homebrew/opt/openblas")
    if(APPLE AND EXISTS "${_openblas_prefix}/lib/libopenblas.dylib")
        list(APPEND GCTA_BLAS_LIBS "${_openblas_prefix}/lib/libopenblas.dylib")
        list(APPEND GCTA_BLAS_INCLUDE_DIRS "${_openblas_prefix}/include")
    else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        list(APPEND GCTA_BLAS_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        if(OPENBLAS_INCLUDE_DIR)
            list(APPEND GCTA_BLAS_INCLUDE_DIRS ${OPENBLAS_INCLUDE_DIR})
        else()
            find_path(LAPACKE_INCLUDE_DIR lapacke.h
                HINTS
                    /usr/include/openblas
                    /usr/local/include/openblas
                    /usr/include/x86_64-linux-gnu/openblas-pthread
                    /usr/include
                    /usr/local/include)
            if(LAPACKE_INCLUDE_DIR)
                list(APPEND GCTA_BLAS_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIR})
                message(STATUS "Found lapacke.h in: ${LAPACKE_INCLUDE_DIR}")
            else()
                message(WARNING "lapacke.h not found. Set OPENBLAS_INCLUDE_DIR to the directory containing lapacke.h")
            endif()
        endif()
    endif()
    list(APPEND GCTA_BLAS_DEFINES GCTA_USE_OPENBLAS EIGEN_USE_LAPACKE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
    set(GCTA_BLAS_PROVIDER "OpenBLAS" CACHE INTERNAL "")

elseif(_GCTA_BLAS STREQUAL "MKL")
    if(LINUX)
        if(NOT OMP_LIBRARY)
            set(OMP_LIBRARY "/usr/lib/libgomp.so")
        endif()
        if(NOT MKL_DIR)
            set(MKL_DIR "/opt/intel/mkl/lib/intel64/cmake/mkl")
        endif()
        set(MKL_INTERFACE "lp64")
        set(MKL_THREADING "gnu_thread")
    endif()
    if(WIN32)
        if(NOT MKL_DIR)
            set(MKL_DIR "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
        endif()
        set(MKL_INTERFACE "lp64")
    endif()
    find_package(MKL REQUIRED)
    list(APPEND GCTA_BLAS_LIBS MKL::MKL)
    list(APPEND GCTA_BLAS_INCLUDE_DIRS ${MKL_INCLUDE})
    list(APPEND GCTA_BLAS_DEFINES GCTA_USE_MKL EIGEN_USE_MKL EIGEN_USE_LAPACKE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
    set(GCTA_BLAS_PROVIDER "MKL" CACHE INTERNAL "")

else()
    message(FATAL_ERROR "Unknown GCTA_BLAS_BACKEND: '${GCTA_BLAS_BACKEND}'. Use Auto, MKL, OpenBLAS, or Accelerate.")
endif()

message(STATUS "BLAS/LAPACK provider: ${GCTA_BLAS_PROVIDER}")

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/include)
set(MAIN_SOURCE "${CMAKE_SOURCE_DIR}/src/main.cpp")

file(GLOB SRCS CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")

# acat.c is pure C — keep it in its own small library so the C++ PCH
# applied to gcta_core does not bleed into C compilation units.
add_library(acat_c OBJECT "${CMAKE_SOURCE_DIR}/src/acat.c")
target_include_directories(acat_c PUBLIC ${COMMON_INCLUDES})

add_library(gcta_core OBJECT ${SRCS})

target_include_directories(gcta_core PUBLIC
    ${COMMON_INCLUDES}
    ${GCTA_BLAS_INCLUDE_DIRS}
    "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include"
    "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/simde")

target_link_libraries(gcta_core PRIVATE
    Boost::algorithm Boost::math Boost::crc Boost::iostreams
    Eigen3::Eigen
    Spectra
    sqlite3
    ${GCTA_ZSTD_LIBS})

target_precompile_headers(gcta_core PRIVATE
    <Eigen/Dense>
    <boost/algorithm/string.hpp>
    <boost/math/distributions.hpp>)

set(GCTA_LIBS gcta_core acat_c)

add_executable(gcta64 ${MAIN_SOURCE})

add_subdirectory(main)
target_include_directories(mainV1 PUBLIC ${COMMON_INCLUDES} ${GCTA_BLAS_INCLUDE_DIRS})
target_link_libraries(mainV1 PRIVATE Eigen3::Eigen Boost::math Boost::iostreams Spectra)

set(GCTA_TARGETS gcta64 mainV1 ${GCTA_LIBS})

set(GCTA_IPO_ENABLED OFF)
if(GCTA_ENABLE_IPO AND CMAKE_BUILD_TYPE MATCHES "^[Rr]elease$")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT GCTA_IPO_SUPPORTED OUTPUT GCTA_IPO_OUTPUT)
    if(GCTA_IPO_SUPPORTED)
        set(GCTA_IPO_ENABLED ON)
    else()
        message(STATUS "IPO/LTO not supported: ${GCTA_IPO_OUTPUT}")
    endif()
endif()

foreach(target IN LISTS GCTA_TARGETS)
    set_target_properties(${target} PROPERTIES
        UNITY_BUILD ${GCTA_UNITY_BUILD}
        INTERPROCEDURAL_OPTIMIZATION ${GCTA_IPO_ENABLED})
    target_include_directories(${target} PRIVATE
        ${GCTA_ZSTD_INCLUDE_DIRS}
        ${GCTA_ZLIB_INCLUDE_DIRS}
        ${GCTA_BLAS_INCLUDE_DIRS}
        $<$<NOT:$<BOOL:${WIN32}>>:${GSL_INCLUDE_DIRS}>)
    target_compile_definitions(${target} PRIVATE
        ${GCTA_BLAS_DEFINES}
        GCTA_BLAS_BACKEND_STR="${GCTA_BLAS_PROVIDER}")
    target_compile_options(${target} PRIVATE
        $<$<BOOL:${APPLE}>:-Wno-nullability-completeness>)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
endforeach()

# Explicit LTO flags per compiler (supplements INTERPROCEDURAL_OPTIMIZATION)
if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease$")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # Thin LTO: most of the benefit of full LTO at a fraction of link time
        target_compile_options(gcta64 PRIVATE -flto=thin)
        target_link_options(gcta64 PRIVATE -flto=thin)
    elseif(CMAKE_COMPILER_IS_GNUCC)
        # -flto=auto uses the same job count as make, no need to query core count
        target_compile_options(gcta64 PRIVATE -flto=auto)
        target_link_options(gcta64 PRIVATE -flto=auto)
    endif()
endif()

add_subdirectory(third_party)

target_include_directories(Pgenlib PRIVATE ${GCTA_ZLIB_INCLUDE_DIRS})

target_link_libraries(gcta64 PRIVATE mainV1 gcta_core acat_c Pgenlib sqlite3 ${GCTA_ZSTD_LIBS} Threads::Threads Eigen3::Eigen)

# Link BLAS before GSL so OpenBLAS/Accelerate/MKL cblas symbols take precedence
# over libgslcblas in dyld (macOS) and ld.so (Linux) symbol resolution order
target_link_libraries(gcta64 PRIVATE ${GCTA_BLAS_LIBS})

if (NOT WIN32)
    # Link only GSL (not gslcblas) when a BLAS backend is available,
    # so GSL uses our BLAS for cblas symbols rather than its own scalar cblas.
    # On Linux wrap in --as-needed/--no-as-needed to prevent GSL from pulling
    # in its own transitive libopenblas which would override our BLAS backend.
    if(_GCTA_BLAS STREQUAL "None")
        target_link_libraries(gcta64 PRIVATE GSL::gsl GSL::gslcblas)
    elseif(UNIX AND NOT APPLE)
        target_link_libraries(gcta64 PRIVATE -Wl,--as-needed GSL::gsl -Wl,--no-as-needed)
    else()
        target_link_libraries(gcta64 PRIVATE GSL::gsl)
    endif()
endif()

# On Linux, re-link BLAS after GSL to ensure it appears last in NEEDED and
# its cblas symbols take precedence via ld.so's first-match ordering.
if(UNIX AND NOT APPLE AND NOT _GCTA_BLAS STREQUAL "None")
    target_link_libraries(gcta64 PRIVATE ${GCTA_BLAS_LIBS})
endif()

foreach(target IN LISTS GCTA_TARGETS)
    if(NOT target STREQUAL "gcta64")
        target_link_libraries(${target} PRIVATE ${GCTA_BLAS_LIBS})
    endif()
endforeach()

# TODO: testing disabled pending fixes
# enable_testing()
# add_subdirectory(test)


# List all targets (verbose, only if requested)
if(CMAKE_VERBOSE_MAKEFILE OR GCTA_LIST_TARGETS)
    foreach(prop IN ITEMS IMPORTED_TARGETS BUILDSYSTEM_TARGETS)
        get_property(_targets DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY ${prop})
        foreach(t IN LISTS _targets)
            message(STATUS "${prop}: ${t}")
        endforeach()
    endforeach()
endif()

# Installation
if(LINUX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
elseif(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
endif()

if(LINUX OR APPLE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)
    install(TARGETS gcta64 mainV1
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
endif()

if(LINUX)
    install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.desktop
        DESTINATION share/applications)
    install(FILES
        ${CMAKE_SOURCE_DIR}/resources/gcta64.png
        DESTINATION share/icons/hicolor/256x256/apps)
    install(FILES
        ${CMAKE_SOURCE_DIR}/resources/gcta64.svg
        DESTINATION share/icons/hicolor/scalable/apps)
elseif(WIN32)
    install(TARGETS gcta64 mainV1 DESTINATION bin)

    file(GLOB RUNTIME_DLLS "${CMAKE_BINARY_DIR}/*.dll")
    install(FILES
        ${RUNTIME_DLLS}
        ${zlib_BINARY_DIR}/zlib.dll
        # Intel MKL
        ${MKL_ROOT}/bin/mkl_intel_thread.2.dll
        ${MKL_ROOT}/bin/mkl_avx2.2.dll
        ${MKL_ROOT}/bin/mkl_core.2.dll
        ${MKL_ROOT}/../../compiler/latest/bin/libiomp5md.dll
        DESTINATION bin)
endif()

# Build configuration summary
message(STATUS "")
message(STATUS "==========================================")
message(STATUS "GCTA2 Build Configuration Summary")
message(STATUS "==========================================")
message(STATUS "  Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "  BLAS/LAPACK:       ${GCTA_BLAS_PROVIDER}")
message(STATUS "  Compiler:          ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Install prefix:    ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  OpenMP:            ${OpenMP_CXX_VERSION}")
if(NOT WIN32)
    message(STATUS "  GSL:               ${GSL_VERSION}")
endif()
message(STATUS "==========================================")
message(STATUS "")