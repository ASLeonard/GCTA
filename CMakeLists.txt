cmake_minimum_required(VERSION 3.28...4.3)

project(GCTA2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
set(Accelerate_DIR "" CACHE PATH "Exact path to the Accelerate framework directory, for example /System/Library/Frameworks/Accelerate.framework")
set(OpenBLAS_DIR "" CACHE PATH "Exact OpenBLAS CMake config directory, for example /opt/homebrew/opt/openblas/lib/cmake/openblas")
set(MKL_DIR "" CACHE PATH "Exact MKL CMake config directory containing MKLConfig.cmake")
set(Boost_DIR "" CACHE PATH "Exact Boost CMake config directory containing BoostConfig.cmake. If set, Boost will be used from this location instead of fetching it.")
option(GCTA_LIST_TARGETS "List imported and local targets during configure" OFF)
option(GCTA_UNITY_BUILD "Enable unity builds for faster compilation" ON)
option(GCTA_ENABLE_IPO "Enable IPO/LTO in Release builds" ON)
option(GCTA_VERBOSE_BLAS "Print verbose BLAS detection information during configure" ON)

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
if(CMAKE_VERBOSE_MAKEFILE)
    set(FETCHCONTENT_QUIET OFF)
else()
    set(FETCHCONTENT_QUIET ON)
endif()

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

FetchContent_Declare(
    stdexec
    URL https://github.com/NVIDIA/stdexec/archive/refs/tags/nvhpc-25.09.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL)

# zlib must be fetched before Boost so Boost iostreams can find it at configure time.
# Kept outside block() so zlib_POPULATED and zlib_SOURCE_DIR are visible in parent scope.
set(ZLIB_BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(zlib)

# Only use Boost when explicitly provided via Boost_DIR; otherwise fall back to
# fetching Boost with FetchContent. Do NOT probe for a system Boost.

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
    set(ZSTD_BUILD_TESTS    OFF)
    set(ZSTD_BUILD_SHARED   ON)
    set(SPECTRA_BUILD_TESTS OFF)
    set(BUILD_TESTING       OFF)
    set(BOOST_IOSTREAMS_ENABLE_ZLIB ON)
    if(zlib_POPULATED)
        set(ZLIB_ROOT        "${zlib_BINARY_DIR}"                      CACHE PATH "" FORCE)
        set(ZLIB_INCLUDE_DIR "${zlib_SOURCE_DIR};${zlib_BINARY_DIR}"   CACHE PATH "" FORCE)
    endif()

    if(Boost_DIR)
        FetchContent_MakeAvailable(zstd)
    else()
        # These two vars are only meaningful to Boost's own CMake build system
        set(BOOST_ENABLE_CMAKE ON)
        set(BOOST_INCLUDE_LIBRARIES "algorithm;math;crc;iostreams")
        FetchContent_MakeAvailable(zstd Boost)
    endif()
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
FetchContent_MakeAvailable(stdexec)

# If the user supplied a Boost_DIR, prefer that over fetching Boost.
if(Boost_DIR)
    message(STATUS "Using external Boost from ${Boost_DIR}")
    find_package(Boost REQUIRED COMPONENTS iostreams
        CONFIG
        PATHS "${Boost_DIR}"
        NO_DEFAULT_PATH)
    find_package(ZLIB REQUIRED)
endif()

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

# Find OpenMP
if(APPLE)
    # Use manual flags for OpenMP with Apple Clang
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_INCLUDE_DIR "/opt/homebrew/opt/libomp/include")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()

if(WIN32)
    set(OpenMP_omp_LIBRARY "${LLVM_ROOT}/lib/libomp.lib")
    set(OpenMP_CXX_FLAGS "-fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_INCLUDE_DIR "${LLVM_ROOT}/lib/clang/20/include")
endif()

find_package(OpenMP REQUIRED)

set(GCTA_BLAS_LIBS "")
set(GCTA_BLAS_INCLUDE_DIRS "")
set(GCTA_BLAS_DEFINES "")

# Resolve "Auto" to a concrete backend
if(GCTA_BLAS_BACKEND STREQUAL "Auto")
    if(APPLE)
        set(_GCTA_BLAS "Accelerate")
    else()
        set(_GCTA_OPENBLAS_CONFIG_FOUND FALSE)
        foreach(_openblas_hint IN LISTS _GCTA_OPENBLAS_HINTS)
            if(EXISTS "${_openblas_hint}/OpenBLASConfig.cmake" OR EXISTS "${_openblas_hint}/openblas-config.cmake")
                set(_GCTA_OPENBLAS_CONFIG_FOUND TRUE)
                break()
            endif()
        endforeach()
        if(_GCTA_OPENBLAS_CONFIG_FOUND)
            set(_GCTA_BLAS "OpenBLAS")
            message(STATUS "Auto-detected OpenBLAS from CONDA_PREFIX")
        else()
            set(_GCTA_BLAS "MKL")
        endif()
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
    if(NOT Accelerate_DIR)
        message(FATAL_ERROR "Accelerate_DIR is required. Set it to the framework directory, for example /System/Library/Frameworks/Accelerate.framework")
    endif()
    list(APPEND GCTA_BLAS_LIBS "${Accelerate_DIR}")
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
    if(NOT OpenBLAS_DIR)
        message(FATAL_ERROR "OpenBLAS_DIR is required. Set it to the OpenBLAS CMake config directory, for example /opt/homebrew/opt/openblas/lib/cmake/openblas")
    endif()
    find_package(OpenBLAS REQUIRED CONFIG PATHS "${OpenBLAS_DIR}" NO_DEFAULT_PATH)

    # Some conda OpenBLAS CMake configs omit OpenBLAS_INCLUDE_DIRS from the
    # exported targets, so cblas.h / lapacke.h won't be on the include path
    # and EIGEN_USE_LAPACKE silently becomes a no-op. Derive the prefix from
    # OpenBLAS_DIR (<prefix>/lib/cmake/openblas) and add it explicitly.
    if(NOT OpenBLAS_INCLUDE_DIRS)
        get_filename_component(_openblas_prefix "${OpenBLAS_DIR}/../../.." ABSOLUTE)
        set(OpenBLAS_INCLUDE_DIRS "${_openblas_prefix}/include")
        message(STATUS "Derived OpenBLAS_INCLUDE_DIRS: ${OpenBLAS_INCLUDE_DIRS}")
        unset(_openblas_prefix)
    endif()
    list(APPEND GCTA_BLAS_INCLUDE_DIRS "${OpenBLAS_INCLUDE_DIRS}")

    if(NOT TARGET OpenBLAS::OpenBLAS)
        if(APPLE)
            add_library(OpenBLAS::OpenBLAS SHARED IMPORTED)
            set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
                IMPORTED_LOCATION "${OpenBLAS_DIR}/../../libopenblas.dylib"
                INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIRS}")
        else()
            # OpenBLAS_DIR = <prefix>/lib/cmake/openblas, so the lib dir is <prefix>/lib
            find_library(OpenBLAS_SO openblas
                HINTS "${OpenBLAS_DIR}/../.."
                NO_DEFAULT_PATH)
            if(OpenBLAS_SO)
                add_library(OpenBLAS::OpenBLAS SHARED IMPORTED)
                set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_SO}"
                    INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIRS}")
            else()
                add_library(OpenBLAS::OpenBLAS STATIC IMPORTED)
                set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
                    IMPORTED_LOCATION "${OpenBLAS_DIR}/../../libopenblas.a"
                    INTERFACE_INCLUDE_DIRECTORIES "${OpenBLAS_INCLUDE_DIRS}")
                find_library(GFORTRAN_LIBRARY gfortran
                    HINTS "${OpenBLAS_DIR}/../../lib")
                if(GFORTRAN_LIBRARY)
                    list(APPEND GCTA_BLAS_LIBS ${GFORTRAN_LIBRARY})
                else()
                    message(FATAL_ERROR "libgfortran not found, required for static OpenBLAS linking")
                endif()
            endif()
        endif()
    endif()
    if(LINUX)
        get_target_property(_openblas_loc OpenBLAS::OpenBLAS IMPORTED_LOCATION)
        if(_openblas_loc)
            get_filename_component(_openblas_libdir "${_openblas_loc}" DIRECTORY)
            list(APPEND CMAKE_BUILD_RPATH "${_openblas_libdir}")
            list(APPEND CMAKE_INSTALL_RPATH "${_openblas_libdir}")
        endif()
    endif()
    list(APPEND GCTA_BLAS_LIBS OpenBLAS::OpenBLAS)
    list(APPEND GCTA_BLAS_DEFINES GCTA_USE_OPENBLAS EIGEN_USE_LAPACKE EIGEN_USE_BLAS EIGEN_USE_LAPACK)
    set(GCTA_BLAS_PROVIDER "OpenBLAS" CACHE INTERNAL "")

elseif(_GCTA_BLAS STREQUAL "MKL")
    if(NOT MKL_DIR)
        message(FATAL_ERROR "MKL_DIR is required. Set it to the MKL CMake config directory containing MKLConfig.cmake")
    endif()
    message(STATUS "Using exact MKL config directory: ${MKL_DIR}")
    # Intel's MKLConfig.cmake defaults the threading layer to "intel_thread"
    # (libiomp5). On Linux this loads a second OpenMP runtime alongside the
    # libgomp already pulled in by find_package(OpenMP), which serialises all
    # parallel regions and causes extreme slowdown. Force "gnu_thread" on Linux
    # so MKL uses the same libgomp as the rest of the build. Let the user
    # override via -DMKL_THREADING=... if they know what they are doing.
    if(NOT DEFINED MKL_THREADING)
        if(LINUX)
            set(MKL_THREADING "gnu_thread")
        else()
            set(MKL_THREADING "intel_thread")
        endif()
    endif()
    if(NOT DEFINED MKL_LINK)
        set(MKL_LINK "dynamic")
    endif()
    if(NOT DEFINED MKL_INTERFACE)
        set(MKL_INTERFACE "lp64")
    endif()
    find_package(MKL REQUIRED CONFIG PATHS "${MKL_DIR}" NO_DEFAULT_PATH)
    # Add the MKL library directory to RPATH on Linux so the binary finds MKL
    # at runtime without relying on LD_LIBRARY_PATH (important in conda envs).
    # MKL_ROOT is set by MKLConfig.cmake; libs live in <root>/lib or
    # <root>/lib/intel64 depending on the installation type.
    if(LINUX AND DEFINED MKL_ROOT)
        foreach(_mkl_lib_subdir IN ITEMS "lib" "lib/intel64")
            set(_mkl_lib_path "${MKL_ROOT}/${_mkl_lib_subdir}")
            if(EXISTS "${_mkl_lib_path}")
                list(APPEND CMAKE_BUILD_RPATH "${_mkl_lib_path}")
                list(APPEND CMAKE_INSTALL_RPATH "${_mkl_lib_path}")
            endif()
        endforeach()
        unset(_mkl_lib_path)
        unset(_mkl_lib_subdir)
    endif()
    list(APPEND GCTA_BLAS_LIBS MKL::MKL)
    list(APPEND GCTA_BLAS_INCLUDE_DIRS ${MKL_INCLUDE})
    list(APPEND GCTA_BLAS_DEFINES GCTA_USE_MKL EIGEN_USE_MKL_ALL)
    set(GCTA_BLAS_PROVIDER "MKL" CACHE INTERNAL "")

else()
    message(FATAL_ERROR "Unknown GCTA_BLAS_BACKEND: '${GCTA_BLAS_BACKEND}'. Use Auto, MKL, OpenBLAS, or Accelerate.")
endif()

message(STATUS "BLAS/LAPACK provider: ${GCTA_BLAS_PROVIDER}")

# Verbose BLAS logging: print detection details when requested
if(CMAKE_VERBOSE_MAKEFILE OR GCTA_VERBOSE_BLAS)
    message(STATUS "---- BLAS detection details ----")
    message(STATUS "  GCTA_BLAS_BACKEND: ${GCTA_BLAS_BACKEND}")
    message(STATUS "  Resolved BLAS backend: ${_GCTA_BLAS}")
    message(STATUS "  OpenBLAS_DIR: ${OpenBLAS_DIR}")
    message(STATUS "  MKL_DIR: ${MKL_DIR}")
    message(STATUS "  MKL_THREADING: ${MKL_THREADING}")
    message(STATUS "  MKL_LINK: ${MKL_LINK}")
    message(STATUS "  MKL_INTERFACE: ${MKL_INTERFACE}")
    message(STATUS "  Accelerate_DIR: ${Accelerate_DIR}")
    message(STATUS "  GCTA_BLAS_INCLUDE_DIRS: ${GCTA_BLAS_INCLUDE_DIRS}")
    message(STATUS "  GCTA_BLAS_DEFINES: ${GCTA_BLAS_DEFINES}")
    message(STATUS "  GCTA_BLAS_LIBS entries: ${GCTA_BLAS_LIBS}")
    foreach(_lib IN LISTS GCTA_BLAS_LIBS)
        if(TARGET ${_lib})
            get_target_property(_loc ${_lib} IMPORTED_LOCATION)
            get_target_property(_incs ${_lib} INTERFACE_INCLUDE_DIRECTORIES)
            message(STATUS "  Target ${_lib}: IMPORTED_LOCATION=${_loc}")
            message(STATUS "  Target ${_lib}: INTERFACE_INCLUDE_DIRECTORIES=${_incs}")
            get_target_property(_type ${_lib} TYPE)
            message(STATUS "  Target ${_lib}: TYPE=${_type}")
        else()
            message(STATUS "  Lib entry ${_lib} (not a CMake target): value='${_lib}'")
        endif()
    endforeach()
    if(TARGET gcta64)
        get_target_property(_gcta_links gcta64 LINK_LIBRARIES)
        message(STATUS "  gcta64 LINK_LIBRARIES: ${_gcta_links}")
        get_target_property(_gcta_linkdirs gcta64 LINK_DIRECTORIES)
        message(STATUS "  gcta64 LINK_DIRECTORIES: ${_gcta_linkdirs}")
    endif()
    message(STATUS "---- end BLAS detection details ----")
endif()

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/include)
set(MAIN_SOURCE "${CMAKE_SOURCE_DIR}/src/main.cpp")

set(SRCS
    ${CMAKE_SOURCE_DIR}/src/Covar.cpp
    ${CMAKE_SOURCE_DIR}/src/FastFAM.cpp
    ${CMAKE_SOURCE_DIR}/src/GRM.cpp
    ${CMAKE_SOURCE_DIR}/src/Geno.cpp
    ${CMAKE_SOURCE_DIR}/src/BedBackend.cpp
    ${CMAKE_SOURCE_DIR}/src/BgenBackend.cpp
    ${CMAKE_SOURCE_DIR}/src/PgenBackend.cpp
    ${CMAKE_SOURCE_DIR}/src/LD.cpp
    ${CMAKE_SOURCE_DIR}/src/Logger.cpp
    ${CMAKE_SOURCE_DIR}/src/Marker.cpp
    ${CMAKE_SOURCE_DIR}/src/OptionIO.cpp
    ${CMAKE_SOURCE_DIR}/src/Pheno.cpp
    ${CMAKE_SOURCE_DIR}/src/StatLib.cpp
    ${CMAKE_SOURCE_DIR}/src/ThreadPool.cpp
    ${CMAKE_SOURCE_DIR}/src/mem.cpp
    ${CMAKE_SOURCE_DIR}/src/tables.cpp
    ${CMAKE_SOURCE_DIR}/src/utils.cpp)

add_executable(gcta64 ${MAIN_SOURCE})

target_sources(gcta64 PRIVATE ${SRCS})

target_include_directories(gcta64 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${COMMON_INCLUDES}
    ${GCTA_BLAS_INCLUDE_DIRS}
    "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include"
    "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/simde")

set(_GCTA_BOOST_LINKS "")
if(TARGET Boost::iostreams)
    list(APPEND _GCTA_BOOST_LINKS Boost::iostreams)
endif()
if(TARGET Boost::headers)
    list(APPEND _GCTA_BOOST_LINKS Boost::headers)
endif()
# Fallback for module-mode find_package (old FindBoost.cmake)
if(DEFINED Boost_LIBRARIES AND NOT "${Boost_LIBRARIES}" STREQUAL "")
    list(APPEND _GCTA_BOOST_LINKS ${Boost_LIBRARIES})
endif()
if(DEFINED Boost_INCLUDE_DIRS)
    target_include_directories(gcta64 PRIVATE ${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(gcta64 PRIVATE
    ${_GCTA_BOOST_LINKS}
    $<$<BOOL:${Boost_DIR}>:ZLIB::ZLIB>   # system Boost::iostreams needs system zlib
    Eigen3::Eigen
    Spectra
    STDEXEC::stdexec
    ${GCTA_ZSTD_LIBS})

target_precompile_headers(gcta64 PRIVATE
    <Eigen/Dense>
    <boost/algorithm/string.hpp>
    <boost/math/distributions.hpp>)

set(GCTA_TARGETS gcta64)

target_sources(gcta64 PRIVATE
    main/CommFunc.cpp
    main/StatFunc.cpp
    main/StrFunc.cpp
    main/bivar_reml.cpp
    main/data.cpp
    main/edata.cpp
    main/ejma.cpp
    main/est_hsq.cpp
    main/gbat.cpp
    main/grm.cpp
    main/gsmr.cpp
    main/gwas_simu.cpp
    main/joint_meta.cpp
    main/ld.cpp
    main/LinAlg.cpp
    main/mbat.cpp
    main/mlm_assoc.cpp
    main/mtcojo.cpp
    main/option.cpp
    main/pc_adjust.cpp
    main/popu_genet.cpp
    main/raw_geno.cpp
    main/reml_within_family.cpp
    main/sbat.cpp)

set(GCTA_TARGETS gcta64)

set(GCTA_IPO_ENABLED OFF)
if(GCTA_ENABLE_IPO AND CMAKE_BUILD_TYPE MATCHES "^[Rr]elease$")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT GCTA_IPO_SUPPORTED OUTPUT GCTA_IPO_OUTPUT LANGUAGES CXX)
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
        ${GCTA_BLAS_INCLUDE_DIRS})
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

add_library(Pgenlib
    third_party/Pgenlib/PgenReader.cpp
    third_party/plink-ng/2.0/include/pgenlib_ffi_support.cc
    third_party/plink-ng/2.0/include/pgenlib_misc.cc
    third_party/plink-ng/2.0/include/pgenlib_read.cc
    third_party/plink-ng/2.0/include/pgenlib_write.cc
    third_party/plink-ng/2.0/include/plink2_base.cc
    third_party/plink-ng/2.0/include/plink2_bgzf.cc
    third_party/plink-ng/2.0/include/plink2_bitmap.cc
    third_party/plink-ng/2.0/include/plink2_bits.cc
    third_party/plink-ng/2.0/include/plink2_fmath.cc
    third_party/plink-ng/2.0/include/plink2_htable.cc
    third_party/plink-ng/2.0/include/plink2_memory.cc
    third_party/plink-ng/2.0/include/plink2_stats.cc
    third_party/plink-ng/2.0/include/plink2_string.cc
    third_party/plink-ng/2.0/include/plink2_text.cc
    third_party/plink-ng/2.0/include/plink2_thread.cc
    third_party/plink-ng/2.0/include/plink2_zstfile.cc
    third_party/plink-ng/2.0/include/pvar_ffi_support.cc)
# Unity builds must be disabled for Pgenlib: it is a separate static library
# with its own include paths and must never be merged into a gcta64 unity TU.
# Without this, CMake can concatenate PgenReader.cpp and PgenBackend.cpp into
# the same translation unit, causing PgenReader.cpp to fail to find GCTA
# headers such as GenoBackendFactory.h.
set_target_properties(Pgenlib PROPERTIES UNITY_BUILD OFF)
target_include_directories(Pgenlib PUBLIC "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include")
if(NOT APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    target_compile_options(Pgenlib PRIVATE -mavx2 -mbmi -mbmi2 -mlzcnt)
endif()

target_include_directories(Pgenlib PRIVATE ${GCTA_ZLIB_INCLUDE_DIRS})

target_link_libraries(gcta64 PRIVATE Pgenlib sqlite3 Threads::Threads)

# Link BLAS early so OpenBLAS/Accelerate/MKL cblas symbols take precedence
# over other cblas implementations in dyld (macOS) and ld.so (Linux) symbol resolution order
target_link_libraries(gcta64 PRIVATE ${GCTA_BLAS_LIBS})

# On Linux, re-link BLAS to ensure it appears last in NEEDED and
# its cblas symbols take precedence via ld.so's first-match ordering.
if(UNIX AND NOT APPLE AND NOT _GCTA_BLAS STREQUAL "None")
    target_link_libraries(gcta64 PRIVATE ${GCTA_BLAS_LIBS})
endif()

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
    # PREPEND so $ORIGIN/../lib is first but BLAS-specific paths appended
    # earlier (e.g. the conda OpenBLAS or MKL lib dir) are preserved.
    # Using set() here would overwrite those previously-appended entries.
    list(PREPEND CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
elseif(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
endif()

if(LINUX OR APPLE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)
    install(TARGETS gcta64
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
    install(TARGETS gcta64 DESTINATION bin)

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

message(STATUS "==========================================")
message(STATUS "")

set_property(DIRECTORY PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/main"
    "${CMAKE_SOURCE_DIR}/third_party/Pgenlib")