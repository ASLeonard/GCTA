cmake_minimum_required(VERSION 3.28...4.3)

project(GCTAg LANGUAGES C CXX)

# Prevent CMake from finding stale configs in build-tree or user package registries.
# This is especially important on HPC clusters where previous partial builds can
# leave broken ZLIBConfig.cmake / zstd-config.cmake files in the build directory.
set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY        TRUE)
set(CMAKE_FIND_PACKAGE_NO_SYSTEM_PACKAGE_REGISTRY TRUE)

# Define platform booleans early — used throughout including in compiler detection.
if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Compiling with GCC ${CMAKE_CXX_COMPILER_VERSION}")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13.0)
        message(FATAL_ERROR "GCTA requires GCC >= 13 for reliable C++23 support (found ${CMAKE_CXX_COMPILER_VERSION})")
    endif()
    # For GCC LTO builds, gcc-ar/gcc-ranlib use the LTO plugin to make archiving
    # LTO-aware (enables cross-archive inlining).  CMake auto-sets
    # CMAKE_CXX_COMPILER_AR for the gcta64 executable when IPO is ON, but static
    # libraries (Pgenlib, sqlite3) use CMAKE_AR — so we need to set it here too.
    # Probe next to the compiler binary so a non-default GCC install is handled
    # correctly even if the system PATH points at a different gcc-ar.
    if(LINUX)
        get_filename_component(_gcc_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
        find_program(_gcc_ar     NAMES gcc-ar     HINTS "${_gcc_bin}" NO_DEFAULT_PATH)
        find_program(_gcc_ranlib NAMES gcc-ranlib HINTS "${_gcc_bin}" NO_DEFAULT_PATH)
        find_program(_gcc_nm     NAMES gcc-nm     HINTS "${_gcc_bin}" NO_DEFAULT_PATH)
        if(_gcc_ar)
            set(CMAKE_AR                  "${_gcc_ar}" CACHE FILEPATH "Archiver" FORCE)
            set(CMAKE_CXX_COMPILER_AR     "${_gcc_ar}" CACHE FILEPATH "" FORCE)
            set(CMAKE_C_COMPILER_AR       "${_gcc_ar}" CACHE FILEPATH "" FORCE)
            message(STATUS "GCC: using LTO-aware archiver ${_gcc_ar}")
        else()
            message(WARNING "GCC: gcc-ar not found next to ${CMAKE_CXX_COMPILER} — LTO may be suboptimal for static libs")
        endif()
        if(_gcc_ranlib)
            set(CMAKE_RANLIB              "${_gcc_ranlib}" CACHE FILEPATH "Ranlib" FORCE)
            set(CMAKE_CXX_COMPILER_RANLIB "${_gcc_ranlib}" CACHE FILEPATH "" FORCE)
            set(CMAKE_C_COMPILER_RANLIB   "${_gcc_ranlib}" CACHE FILEPATH "" FORCE)
        endif()
        if(_gcc_nm)
            set(CMAKE_NM                  "${_gcc_nm}" CACHE FILEPATH "NM" FORCE)
            set(CMAKE_CXX_COMPILER_NM     "${_gcc_nm}" CACHE FILEPATH "" FORCE)
            set(CMAKE_C_COMPILER_NM       "${_gcc_nm}" CACHE FILEPATH "" FORCE)
        endif()
        # Pin libgomp to this GCC's copy so HPC LD_LIBRARY_PATH doesn't shadow it.
        # Uses the same -print-file-name probe as the AOCC libomp / toolchain logic.
        execute_process(
            COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libgomp.so.1
            OUTPUT_VARIABLE _gcc_gomp_path
            OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        if(_gcc_gomp_path AND NOT _gcc_gomp_path STREQUAL "libgomp.so.1")
            file(REAL_PATH "${_gcc_gomp_path}" _gcc_gomp_real)
            get_filename_component(_gcc_gomp_dir "${_gcc_gomp_real}" DIRECTORY)
            list(APPEND CMAKE_BUILD_RPATH   "${_gcc_gomp_dir}")
            list(APPEND CMAKE_INSTALL_RPATH "${_gcc_gomp_dir}")
            message(STATUS "GCC: libgomp.so.1 -> ${_gcc_gomp_real}")
        endif()
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Compiling with Clang ${CMAKE_CXX_COMPILER_VERSION}")
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

# Detect AOCC (AMD's Clang-based compiler). It reports CMAKE_CXX_COMPILER_ID
# as "Clang" but has its own libomp and a clang-specific AOCL library variant.
# We probe the --version string since there is no dedicated CMake variable.
set(GCTA_COMPILER_IS_AOCC FALSE)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND LINUX)
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} --version
        OUTPUT_VARIABLE _cxx_ver_str
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(_cxx_ver_str MATCHES "AOCC")
        set(GCTA_COMPILER_IS_AOCC TRUE)
        # Derive the AOCC install root from the compiler path (bin/../)
        get_filename_component(_aocc_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
        get_filename_component(GCTA_AOCC_ROOT "${_aocc_bin}" DIRECTORY CACHE)
        message(STATUS "Detected AOCC compiler (root: ${GCTA_AOCC_ROOT})")

        # AOCC requires its own llvm-ar / llvm-ranlib for static library creation.
        # The system 'ar' cannot handle LLVM bitcode objects produced by AOCC,
        # causing "Error running link command: no such file or directory" at archive time.
        find_program(_aocc_ar     NAMES llvm-ar     HINTS "${GCTA_AOCC_ROOT}/bin" NO_DEFAULT_PATH)
        find_program(_aocc_ranlib NAMES llvm-ranlib HINTS "${GCTA_AOCC_ROOT}/bin" NO_DEFAULT_PATH)
        find_program(_aocc_nm     NAMES llvm-nm     HINTS "${GCTA_AOCC_ROOT}/bin" NO_DEFAULT_PATH)
        if(_aocc_ar)
            set(CMAKE_AR               "${_aocc_ar}" CACHE FILEPATH "Archiver" FORCE)
            set(CMAKE_CXX_COMPILER_AR  "${_aocc_ar}" CACHE FILEPATH "" FORCE)
            set(CMAKE_C_COMPILER_AR    "${_aocc_ar}" CACHE FILEPATH "" FORCE)
            message(STATUS "AOCC: using archiver ${_aocc_ar}")
        else()
            message(WARNING "AOCC: llvm-ar not found under ${GCTA_AOCC_ROOT}/bin — static library builds may fail")
        endif()
        if(_aocc_ranlib)
            set(CMAKE_RANLIB               "${_aocc_ranlib}" CACHE FILEPATH "Ranlib" FORCE)
            set(CMAKE_CXX_COMPILER_RANLIB  "${_aocc_ranlib}" CACHE FILEPATH "" FORCE)
            set(CMAKE_C_COMPILER_RANLIB    "${_aocc_ranlib}" CACHE FILEPATH "" FORCE)
        endif()
        if(_aocc_nm)
            set(CMAKE_NM               "${_aocc_nm}" CACHE FILEPATH "NM" FORCE)
            set(CMAKE_CXX_COMPILER_NM  "${_aocc_nm}" CACHE FILEPATH "" FORCE)
            set(CMAKE_C_COMPILER_NM    "${_aocc_nm}" CACHE FILEPATH "" FORCE)
        endif()
    endif()
endif()

# For upstream Clang (not AOCC) on Linux, detect LLVM tooling next to the
# compiler so ar/nm/ranlib come from the same LLVM install rather than PATH.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND LINUX AND NOT GCTA_COMPILER_IS_AOCC)
    get_filename_component(_clang_bin "${CMAKE_CXX_COMPILER}" DIRECTORY)
    find_program(_clang_ar     NAMES llvm-ar     HINTS "${_clang_bin}" NO_DEFAULT_PATH)
    find_program(_clang_ranlib NAMES llvm-ranlib HINTS "${_clang_bin}" NO_DEFAULT_PATH)
    find_program(_clang_nm     NAMES llvm-nm     HINTS "${_clang_bin}" NO_DEFAULT_PATH)
    if(_clang_ar)
        set(CMAKE_AR               "${_clang_ar}" CACHE FILEPATH "Archiver" FORCE)
        set(CMAKE_CXX_COMPILER_AR  "${_clang_ar}" CACHE FILEPATH "" FORCE)
        set(CMAKE_C_COMPILER_AR    "${_clang_ar}" CACHE FILEPATH "" FORCE)
        message(STATUS "Clang: using archiver ${_clang_ar}")
    endif()
    if(_clang_ranlib)
        set(CMAKE_RANLIB               "${_clang_ranlib}" CACHE FILEPATH "Ranlib" FORCE)
        set(CMAKE_CXX_COMPILER_RANLIB  "${_clang_ranlib}" CACHE FILEPATH "" FORCE)
        set(CMAKE_C_COMPILER_RANLIB    "${_clang_ranlib}" CACHE FILEPATH "" FORCE)
    endif()
    if(_clang_nm)
        set(CMAKE_NM               "${_clang_nm}" CACHE FILEPATH "NM" FORCE)
        set(CMAKE_CXX_COMPILER_NM  "${_clang_nm}" CACHE FILEPATH "" FORCE)
        set(CMAKE_C_COMPILER_NM    "${_clang_nm}" CACHE FILEPATH "" FORCE)
    endif()
endif()

# AOCC is LLVM 17 and cannot parse headers from GCC 15+. If your default GCC is
# too new (check: gcc --version), point AOCC at an older GCC toolchain:
#   -DGCTA_GCC_TOOLCHAIN=/path/to/gcc-13-or-14-root
# The path should be the GCC prefix — the directory containing bin/gcc.
set(GCTA_GCC_TOOLCHAIN "" CACHE PATH
    "GCC installation root for AOCC to use as its C++ header/sysroot provider. \
Required when the system GCC is >= 15 (AOCC/Clang 17 cannot parse GCC 15+ libstdc++ headers). \
Example: /opt/gcc-13  (must contain bin/gcc)")
if(GCTA_COMPILER_IS_AOCC AND GCTA_GCC_TOOLCHAIN)
    if(NOT IS_DIRECTORY "${GCTA_GCC_TOOLCHAIN}")
        message(FATAL_ERROR "GCTA_GCC_TOOLCHAIN does not exist: ${GCTA_GCC_TOOLCHAIN}")
    endif()
    message(STATUS "AOCC: using GCC toolchain at ${GCTA_GCC_TOOLCHAIN}")
    add_compile_options(--gcc-toolchain=${GCTA_GCC_TOOLCHAIN})
    add_link_options(--gcc-toolchain=${GCTA_GCC_TOOLCHAIN})

    # Use the compiler's own -print-file-name to locate libstdc++.so.6 and
    # libgcc_s.so.1 within the toolchain.  This is more reliable than directory
    # globs: the compiler applies the same search logic (multilib, sysroot,
    # --gcc-toolchain) it uses at link time, so we get the exact file it will
    # link against.  Spack and other toolchain managers emit paths with ../../../../
    # traversals (e.g. lib/gcc/<triple>/<ver>/../../../../lib64/); file(REAL_PATH)
    # normalises these to canonical absolute paths before adding to RPATH.
    foreach(_tc_lib IN ITEMS libstdc++.so.6 libgcc_s.so.1)
        execute_process(
            COMMAND ${CMAKE_CXX_COMPILER}
                    --gcc-toolchain=${GCTA_GCC_TOOLCHAIN}
                    -print-file-name=${_tc_lib}
            OUTPUT_VARIABLE _tc_lib_path
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        # The compiler echoes back the bare filename when it cannot find the file.
        if(_tc_lib_path AND NOT _tc_lib_path STREQUAL "${_tc_lib}")
            file(REAL_PATH "${_tc_lib_path}" _tc_lib_real)
            get_filename_component(_tc_lib_dir "${_tc_lib_real}" DIRECTORY)
            list(APPEND CMAKE_BUILD_RPATH   "${_tc_lib_dir}")
            list(APPEND CMAKE_INSTALL_RPATH "${_tc_lib_dir}")
            message(STATUS "AOCC: ${_tc_lib} -> ${_tc_lib_real}")
        else()
            message(WARNING
                "AOCC: -print-file-name=${_tc_lib} did not resolve to a path "
                "in ${GCTA_GCC_TOOLCHAIN} — runtime may load the system library. "
                "Use -DGCTA_STATIC_STDCXX=ON to avoid the runtime dependency.")
        endif()
    endforeach()
    message(STATUS "AOCC: added GCC toolchain lib dirs to rpath")
elseif(GCTA_COMPILER_IS_AOCC AND NOT GCTA_GCC_TOOLCHAIN)
    # Warn if the detected system GCC is too new for AOCC to handle.
    execute_process(
        COMMAND gcc --version
        OUTPUT_VARIABLE _gcc_ver_str
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(_gcc_ver_str MATCHES "([0-9]+)\\.[0-9]+\\.[0-9]+")
        if(CMAKE_MATCH_1 GREATER_EQUAL 15)
            message(WARNING
                "AOCC detected but system GCC ${CMAKE_MATCH_1} is too new — "
                "AOCC (Clang 17) cannot parse GCC 15+ libstdc++ headers.\n"
                "Set -DGCTA_GCC_TOOLCHAIN=/path/to/gcc-13-or-14 to fix compile errors.")
        endif()
    endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "Defaulting to Release build")
endif()

# Release-build arch flags are applied per-target on gcta64 below,
# after add_executable, so they don't bleed into FetchContent dependencies.

set(GCTA_BLAS_BACKEND "Auto" CACHE STRING "BLAS/LAPACK backend: Auto, OpenBLAS, MKL, Accelerate, or AOCL")
set_property(CACHE GCTA_BLAS_BACKEND PROPERTY STRINGS Auto OpenBLAS MKL Accelerate AOCL)
set(GCTA_BLAS_LIBRARY "" CACHE FILEPATH "Primary BLAS library (OpenBLAS/MKL only). E.g. /usr/lib/x86_64-linux-gnu/libopenblas.so or /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so")
set(GCTA_LAPACKE_LIBRARY "" CACHE FILEPATH "Separate LAPACKE library, if needed (e.g. /usr/lib/x86_64-linux-gnu/liblapacke.so for apt OpenBLAS). Leave empty if your BLAS library bundles LAPACKE (MKL, conda OpenBLAS).")
set(GCTA_BLAS_INCLUDE_DIR "" CACHE PATH "Directory containing BLAS/LAPACK headers (OpenBLAS/MKL only). E.g. /usr/include/x86_64-linux-gnu/openblas or /opt/intel/oneapi/mkl/latest/include")
set(Boost_DIR "" CACHE PATH "Exact Boost CMake config directory containing BoostConfig.cmake. If set, Boost will be used from this location instead of fetching it.")
option(GCTA_LIST_TARGETS "List imported and local targets during configure" OFF)
option(GCTA_UNITY_BUILD "Enable unity builds for faster compilation" ON)
option(GCTA_ENABLE_IPO "Enable IPO/LTO in Release builds" ON)
option(GCTA_VERBOSE_BLAS "Print verbose BLAS detection information during configure" OFF)
option(BGEN_SUPPORT "Enable BGEN file format support (requires SQLite)" OFF)
option(GCTA_STATIC_STDCXX
    "Statically link libstdc++ and libgcc into the binary (AOCC/Clang on Linux). \
Avoids runtime dependency on the toolchain's shared libstdc++.so.6 / libgcc_s.so.1. \
Recommended when deploying to nodes that do not have the build toolchain's GCC modules loaded."
    OFF)

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
# Once all tarballs are cached locally, skip re-checking download URLs on
# every re-configure.  This is faster and works on HPC compute nodes that
# have no outbound internet.  Override with -DFETCHCONTENT_UPDATES_DISCONNECTED=OFF.
set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL
    "Skip re-checking FetchContent URLs on re-configure (set OFF to force update check)")
if(CMAKE_VERBOSE_MAKEFILE)
    set(FETCHCONTENT_QUIET OFF)
else()
    set(FETCHCONTENT_QUIET ON)
endif()

FetchContent_Declare(
    zlib
    URL https://github.com/madler/zlib/archive/refs/tags/v1.3.2.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL)
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
    EXCLUDE_FROM_ALL)
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.91.0/boost-1.91.0-cmake.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL)
if(BGEN_SUPPORT)
    FetchContent_Declare(
        sqlite
        URL https://www.sqlite.org/2026/sqlite-amalgamation-3510200.zip
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(
    stdexec
    URL https://github.com/NVIDIA/stdexec/archive/refs/tags/nvhpc-26.05.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    EXCLUDE_FROM_ALL)

# zlib must be fetched before Boost so Boost iostreams can find it at configure time.
# Kept outside block() so zlib_POPULATED and zlib_SOURCE_DIR are visible in parent scope.
set(ZLIB_BUILD_EXAMPLES OFF)
FetchContent_MakeAvailable(zlib)

# Only use Boost when explicitly provided via Boost_DIR; otherwise fall back to
# fetching Boost with FetchContent. Do NOT probe for a system Boost.

# After zlib is available, set hints so Boost iostreams links against it correctly.
# The target-based hints below (ZLIB_LIBRARY = ZLIB::ZLIB) are the authoritative
# version; this block is intentionally removed to avoid a stale path-based override.

block()
    set(ZSTD_BUILD_PROGRAMS OFF)
    set(ZSTD_BUILD_TESTS    OFF)
    set(ZSTD_BUILD_SHARED   ON)
    set(SPECTRA_BUILD_TESTS OFF)
    set(BUILD_TESTING       OFF)
    set(BOOST_IOSTREAMS_ENABLE_ZLIB ON)
    if(zlib_POPULATED AND TARGET ZLIB::ZLIB)
        set(ZLIB_ROOT        "${zlib_BINARY_DIR}"                    CACHE PATH   "" FORCE)
        set(ZLIB_LIBRARY     ZLIB::ZLIB                              CACHE STRING "" FORCE)
        set(ZLIB_INCLUDE_DIRS "${zlib_SOURCE_DIR};${zlib_BINARY_DIR}" CACHE PATH   "" FORCE)
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
if(BGEN_SUPPORT)
    FetchContent_MakeAvailable(sqlite)
endif()
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
# zlib is always built from source by FetchContent; alias whichever static
# target its CMakeLists creates ('zlibstatic' on all platforms, or 'zlib' on UNIX).
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

# Extra hints so Boost's internal find_package(ZLIB) also picks up the fetched zlib.
# Using the target directly avoids the configure-time / build-time ordering problem
# (libz.a does not exist yet when CMake configures Boost, but the target does).
if(zlib_POPULATED AND TARGET ZLIB::ZLIB)
    set(ZLIB_ROOT         "${zlib_BINARY_DIR}"                     CACHE PATH   "" FORCE)
    set(ZLIB_LIBRARY      ZLIB::ZLIB                               CACHE STRING "" FORCE)
    set(ZLIB_INCLUDE_DIRS "${zlib_SOURCE_DIR};${zlib_BINARY_DIR}"  CACHE PATH   "" FORCE)
    set(ZLIB_USE_STATIC_LIBS ON                                    CACHE BOOL   "" FORCE)
endif()

# Resolve zstd link target.
# Check for the actual CMake target rather than relying on zstd_POPULATED,
# which CMake 4.x no longer sets reliably after a block() scope.
# FetchContent's zstd build always creates 'libzstd'; system installs use
# the zstd:: namespace targets.
if(TARGET libzstd)
    # Built from source by FetchContent — target carries its own include dirs
    set(GCTA_ZSTD_INCLUDE_DIRS "")
    set(GCTA_ZSTD_LIBS libzstd)
elseif(TARGET zstd::libzstd_shared)
    set(GCTA_ZSTD_INCLUDE_DIRS "")
    set(GCTA_ZSTD_LIBS zstd::libzstd_shared)
elseif(TARGET zstd::libzstd_static)
    set(GCTA_ZSTD_INCLUDE_DIRS "")
    set(GCTA_ZSTD_LIBS zstd::libzstd_static)
else()
    message(FATAL_ERROR
        "zstd is not available: FetchContent did not create a libzstd target "
        "and no system zstd:: targets were found.\n"
        "Try: module unload zstd && rm -rf build && cmake -B build ...")
endif()

if(BGEN_SUPPORT)
    add_library(sqlite3 STATIC ${sqlite_SOURCE_DIR}/sqlite3.c)
    target_include_directories(sqlite3 PUBLIC ${sqlite_SOURCE_DIR})
endif()

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

# AOCC ships its own libomp — steer CMake to it before it finds a system libgomp.
# Without this, find_package(OpenMP) may pick up GCC's libgomp instead, which
# causes runtime crashes when the binary was compiled with AOCC.
if(GCTA_COMPILER_IS_AOCC)
    set(OpenMP_CXX_LIB_NAMES "omp")
    # Look for libomp.so under the AOCC root; AOCC 4.x puts it in lib/
    find_library(_aocc_omp_lib
        NAMES omp
        HINTS "${GCTA_AOCC_ROOT}/lib"
        NO_DEFAULT_PATH)
    if(_aocc_omp_lib)
        set(OpenMP_omp_LIBRARY "${_aocc_omp_lib}")
        message(STATUS "AOCC OpenMP library: ${_aocc_omp_lib}")
        # Ensure the AOCC lib dir is in the rpath so the binary finds libomp at runtime
        get_filename_component(_aocc_omp_dir "${_aocc_omp_lib}" DIRECTORY)
        list(APPEND CMAKE_BUILD_RPATH   "${_aocc_omp_dir}")
        list(APPEND CMAKE_INSTALL_RPATH "${_aocc_omp_dir}")
    else()
        message(WARNING
            "AOCC detected but libomp.so not found under ${GCTA_AOCC_ROOT}/lib.\n"
            "OpenMP detection may fall back to the system libgomp, which can cause "
            "runtime crashes. Set OpenMP_omp_LIBRARY manually if needed.")
    endif()
endif()

find_package(OpenMP REQUIRED)

set(GCTA_BLAS_LIBS "")
set(GCTA_BLAS_INCLUDE_DIRS "")
set(GCTA_BLAS_DEFINES "")

# Resolve "Auto" backend.
# On Linux/Windows the user must supply explicit paths — CMake's FindBLAS linker
# probes are too fragile across distros and package managers to be reliable.
if(GCTA_BLAS_BACKEND STREQUAL "Auto")
    if(APPLE)
        set(_GCTA_BLAS "Accelerate")
    else()
        message(FATAL_ERROR
            "GCTA_BLAS_BACKEND=Auto is not supported on Linux/Windows.\n"
            "Set -DGCTA_BLAS_BACKEND=OpenBLAS or -DGCTA_BLAS_BACKEND=MKL and supply:\n"
            "  -DGCTA_BLAS_LIBRARY=<full path to .so/.a>\n"
            "  -DGCTA_BLAS_INCLUDE_DIR=<directory containing cblas.h / mkl.h>\n"
            "  -DGCTA_LAPACKE_LIBRARY=<separate lapacke .so if needed, else omit>\n"
            "OpenBLAS example (Ubuntu x86-64):\n"
            "  -DGCTA_BLAS_LIBRARY=/usr/lib/x86_64-linux-gnu/libopenblas.so\n"
            "  -DGCTA_LAPACKE_LIBRARY=/usr/lib/x86_64-linux-gnu/liblapacke.so\n"
            "  -DGCTA_BLAS_INCLUDE_DIR=/usr/include/x86_64-linux-gnu/openblas\n"
            "MKL example (set MKL_THREADING_LAYER=GNU with GCC, or MKL_THREADING_LAYER=IOMP5 with AOCC):\n"
            "  -DGCTA_BLAS_LIBRARY=/opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so\n"
            "  -DGCTA_BLAS_INCLUDE_DIR=/opt/intel/oneapi/mkl/latest/include")
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
if(APPLE AND _GCTA_BLAS STREQUAL "AOCL")
    message(FATAL_ERROR "AOCL is not supported on macOS. Use Accelerate or OpenBLAS.")
endif()
if(NOT APPLE AND _GCTA_BLAS STREQUAL "Accelerate")
    message(FATAL_ERROR "Accelerate is only available on macOS. Use MKL, OpenBLAS, or AOCL.")
endif()

if(_GCTA_BLAS STREQUAL "Accelerate")
    # Accelerate is a macOS system framework — no user-supplied paths needed.
    set(BLA_VENDOR Apple)
    find_package(BLAS   REQUIRED)
    find_package(LAPACK REQUIRED)

    if(NOT CMAKE_OSX_SYSROOT)
        execute_process(
            COMMAND xcrun --sdk macosx --show-sdk-path
            OUTPUT_VARIABLE CMAKE_OSX_SYSROOT
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    find_path(GCTA_ACCELERATE_INCLUDE_DIR
        NAMES Accelerate.h cblas.h
        HINTS "${CMAKE_OSX_SYSROOT}/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers"
              "${CMAKE_OSX_SYSROOT}/System/Library/Frameworks/vecLib.framework/Headers"
        NO_DEFAULT_PATH
        DOC "vecLib / Accelerate headers directory")
    if(NOT GCTA_ACCELERATE_INCLUDE_DIR)
        message(FATAL_ERROR "Accelerate/vecLib headers not found inside SDK: ${CMAKE_OSX_SYSROOT}")
    endif()
    message(STATUS "Accelerate include dir: ${GCTA_ACCELERATE_INCLUDE_DIR}")

    list(APPEND GCTA_BLAS_INCLUDE_DIRS "${GCTA_ACCELERATE_INCLUDE_DIR}")
    list(APPEND GCTA_BLAS_LIBS         ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    list(APPEND GCTA_BLAS_DEFINES      GCTA_USE_ACCELERATE ACCELERATE_NEW_LAPACK EIGEN_USE_BLAS)
    set(GCTA_BLAS_PROVIDER "Accelerate" CACHE INTERNAL "")

elseif(_GCTA_BLAS STREQUAL "OpenBLAS" OR _GCTA_BLAS STREQUAL "MKL" OR _GCTA_BLAS STREQUAL "AOCL")
    # OpenBLAS, MKL, and AOCL all require the user to supply explicit paths.
    # This avoids the fragile FindBLAS / BLA_VENDOR linker-probe machinery.
    # For MKL, pass libmkl_rt.so (the single-library interface).
    #   With GCC:  set MKL_THREADING_LAYER=GNU at runtime to use libgomp.
    #   With AOCC: set MKL_THREADING_LAYER=IOMP5 at runtime to use libomp.
    # For AOCL, use the gcc/ library variant with GCC and the clang/ variant with AOCC.
    if(NOT GCTA_BLAS_LIBRARY)
        message(FATAL_ERROR
            "GCTA_BLAS_LIBRARY must be set for the ${_GCTA_BLAS} backend.\n"
            "Examples:\n"
            "  OpenBLAS: -DGCTA_BLAS_LIBRARY=/usr/lib/x86_64-linux-gnu/libopenblas.so\n"
            "            -DGCTA_LAPACKE_LIBRARY=/usr/lib/x86_64-linux-gnu/liblapacke.so\n"
            "  MKL:      -DGCTA_BLAS_LIBRARY=/opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so\n"
            "  AOCL (GCC):  -DGCTA_BLAS_LIBRARY=/opt/AMD/aocl/aocl-linux-gcc-5.0.0/gcc/lib/libblis.so\n"
            "               -DGCTA_LAPACKE_LIBRARY=/opt/AMD/aocl/aocl-linux-gcc-5.0.0/gcc/lib/libflame.so\n"
            "  AOCL (AOCC): -DGCTA_BLAS_LIBRARY=/opt/AMD/aocl/aocl-linux-clang-5.0.0/clang/lib/libblis.so\n"
            "               -DGCTA_LAPACKE_LIBRARY=/opt/AMD/aocl/aocl-linux-clang-5.0.0/clang/lib/libflame.so")
    endif()
    if(NOT GCTA_BLAS_INCLUDE_DIR)
        message(FATAL_ERROR
            "GCTA_BLAS_INCLUDE_DIR must be set for the ${_GCTA_BLAS} backend.\n"
            "Examples:\n"
            "  OpenBLAS:    -DGCTA_BLAS_INCLUDE_DIR=/usr/include/x86_64-linux-gnu/openblas\n"
            "  MKL:         -DGCTA_BLAS_INCLUDE_DIR=/opt/intel/oneapi/mkl/latest/include\n"
            "  AOCL (GCC):  -DGCTA_BLAS_INCLUDE_DIR=/opt/AMD/aocl/aocl-linux-gcc-5.0.0/gcc/include\n"
            "  AOCL (AOCC): -DGCTA_BLAS_INCLUDE_DIR=/opt/AMD/aocl/aocl-linux-clang-5.0.0/clang/include")
    endif()
    if(NOT EXISTS "${GCTA_BLAS_LIBRARY}")
        message(FATAL_ERROR "GCTA_BLAS_LIBRARY does not exist: ${GCTA_BLAS_LIBRARY}")
    endif()
    if(NOT IS_DIRECTORY "${GCTA_BLAS_INCLUDE_DIR}")
        message(FATAL_ERROR "GCTA_BLAS_INCLUDE_DIR does not exist: ${GCTA_BLAS_INCLUDE_DIR}")
    endif()

    message(STATUS "${_GCTA_BLAS} library:     ${GCTA_BLAS_LIBRARY}")
    if(GCTA_LAPACKE_LIBRARY)
        message(STATUS "${_GCTA_BLAS} LAPACKE library: ${GCTA_LAPACKE_LIBRARY}")
    endif()
    message(STATUS "${_GCTA_BLAS} include dir: ${GCTA_BLAS_INCLUDE_DIR}")

    if(LINUX)
        get_filename_component(_blas_lib_dir "${GCTA_BLAS_LIBRARY}" DIRECTORY)
        list(APPEND CMAKE_BUILD_RPATH   "${_blas_lib_dir}")
        list(APPEND CMAKE_INSTALL_RPATH "${_blas_lib_dir}")
        # --rpath-link lets the linker resolve DT_NEEDED entries of your .so deps
        # without those libs needing to be explicit link targets themselves
        add_link_options("LINKER:-rpath-link,${_blas_lib_dir}")
        if(GCTA_LAPACKE_LIBRARY)
            get_filename_component(_lapack_lib_dir "${GCTA_LAPACKE_LIBRARY}" DIRECTORY)
            list(APPEND CMAKE_BUILD_RPATH   "${_lapack_lib_dir}")
            list(APPEND CMAKE_INSTALL_RPATH "${_lapack_lib_dir}")
            add_link_options("LINKER:-rpath-link,${_lapack_lib_dir}")
        endif()
    endif()

    list(APPEND GCTA_BLAS_INCLUDE_DIRS "${GCTA_BLAS_INCLUDE_DIR}")
    list(APPEND GCTA_BLAS_LIBS         "${GCTA_BLAS_LIBRARY}")
    if(GCTA_LAPACKE_LIBRARY)
        if(NOT EXISTS "${GCTA_LAPACKE_LIBRARY}")
            message(FATAL_ERROR "GCTA_LAPACKE_LIBRARY does not exist: ${GCTA_LAPACKE_LIBRARY}")
        endif()
        list(APPEND GCTA_BLAS_LIBS "${GCTA_LAPACKE_LIBRARY}")
    endif()

    if(_GCTA_BLAS STREQUAL "OpenBLAS")
        list(APPEND GCTA_BLAS_DEFINES GCTA_USE_OPENBLAS EIGEN_USE_LAPACKE EIGEN_USE_BLAS)
        set(GCTA_BLAS_PROVIDER "OpenBLAS" CACHE INTERNAL "")
    elseif(_GCTA_BLAS STREQUAL "AOCL")
        list(APPEND GCTA_BLAS_DEFINES GCTA_USE_AOCL EIGEN_USE_LAPACKE EIGEN_USE_BLAS)
        set(GCTA_BLAS_PROVIDER "AOCL" CACHE INTERNAL "")

        # AOCL ships two separate library builds:
        #   aocc/  — linked against libomp  (LLVM OpenMP, AOCC's runtime)
        #   gcc/   — linked against libgomp (GCC's OpenMP runtime)
        # Using the AOCC variant with GCC (or vice-versa) produces link errors:
        #   "libomp.so not found" / "undefined reference to __kmpc_*"
        # Detect the mismatch early via the presence of libomp.so in the lib dir.
        get_filename_component(_aocl_lib_dir "${GCTA_BLAS_LIBRARY}" DIRECTORY)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND EXISTS "${_aocl_lib_dir}/libomp.so")
            message(FATAL_ERROR
                "AOCL compiler/variant mismatch: you are using GCC but the AOCL "
                "library at\n  ${_aocl_lib_dir}\nwas built for AOCC (it contains "
                "libomp.so, LLVM OpenMP).\n"
                "Use the GCC-built AOCL variant instead, e.g.:\n"
                "  -DGCTA_BLAS_LIBRARY=<aocl-root>/gcc/MT/lib/libblis-mt.so\n"
                "  -DGCTA_LAPACKE_LIBRARY=<aocl-root>/gcc/MT/lib/libflame.so\n"
                "  -DGCTA_BLAS_INCLUDE_DIR=<aocl-root>/gcc/MT/include")
        endif()
        if(GCTA_COMPILER_IS_AOCC AND NOT EXISTS "${_aocl_lib_dir}/libomp.so"
                                 AND EXISTS "${_aocl_lib_dir}/libgomp.so")
            message(WARNING
                "AOCL compiler/variant mismatch: you are using AOCC but the AOCL "
                "library at\n  ${_aocl_lib_dir}\nappears to be the GCC-built variant "
                "(no libomp.so found).\nConsider using the AOCC-built variant for "
                "best performance:\n"
                "  -DGCTA_BLAS_LIBRARY=<aocl-root>/aocc/MT/lib/libblis-mt.so")
        endif()

        # AOCL transitive dependencies: libflame needs libaoclutils and libblis-mt.
        # Skip any dep whose path is already in GCTA_BLAS_LIBS (e.g. libblis-mt
        # is the primary BLAS library and must not be linked a second time).
        get_filename_component(_aocl_lib_dir "${GCTA_BLAS_LIBRARY}" DIRECTORY)
        foreach(_aocl_dep IN ITEMS libaoclutils.so libblis-mt.so)
            set(_aocl_dep_path "${_aocl_lib_dir}/${_aocl_dep}")
            if("${_aocl_dep_path}" IN_LIST GCTA_BLAS_LIBS)
                message(STATUS "AOCL: skipping ${_aocl_dep} (already linked as primary library)")
            elseif(EXISTS "${_aocl_dep_path}")
                list(APPEND GCTA_BLAS_LIBS "${_aocl_dep_path}")
                message(STATUS "AOCL: found transitive dep ${_aocl_dep}")
            else()
                message(WARNING "AOCL: transitive dep not found: ${_aocl_dep_path}")
            endif()
        endforeach()
    else()
        list(APPEND GCTA_BLAS_DEFINES GCTA_USE_MKL EIGEN_USE_MKL_ALL)
        set(GCTA_BLAS_PROVIDER "MKL" CACHE INTERNAL "")
    endif()

else()
    message(FATAL_ERROR "Unknown GCTA_BLAS_BACKEND: '${GCTA_BLAS_BACKEND}'. Use Auto, OpenBLAS, MKL, Accelerate, or AOCL.")
endif()

message(STATUS "BLAS/LAPACK provider: ${GCTA_BLAS_PROVIDER}")

# Verbose BLAS logging: print detection details when requested
if(CMAKE_VERBOSE_MAKEFILE OR GCTA_VERBOSE_BLAS)
    message(STATUS "---- BLAS detection details ----")
    message(STATUS "  GCTA_BLAS_BACKEND:    ${GCTA_BLAS_BACKEND}")
    message(STATUS "  Resolved backend:     ${_GCTA_BLAS}")
    message(STATUS "  GCTA_BLAS_LIBRARY:    ${GCTA_BLAS_LIBRARY}")
    message(STATUS "  GCTA_LAPACKE_LIBRARY: ${GCTA_LAPACKE_LIBRARY}")
    message(STATUS "  GCTA_BLAS_INCLUDE_DIR:${GCTA_BLAS_INCLUDE_DIR}")
    message(STATUS "  Defines:              ${GCTA_BLAS_DEFINES}")
    message(STATUS "  Link list:")
    foreach(_lib IN LISTS GCTA_BLAS_LIBS)
        message(STATUS "    ${_lib}")
    endforeach()
    message(STATUS "---- end BLAS detection details ----")
endif()

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/include)
set(MAIN_SOURCE "${CMAKE_SOURCE_DIR}/src/main.cpp")

set(SRCS
    ${CMAKE_SOURCE_DIR}/src/Covar.cpp
    ${CMAKE_SOURCE_DIR}/src/FastFAM.cpp
    ${CMAKE_SOURCE_DIR}/src/MLMA_stream.cpp
    ${CMAKE_SOURCE_DIR}/src/MLMA_loco.cpp
    ${CMAKE_SOURCE_DIR}/src/RemlEngine.cpp
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

# Default OFF; enabled below only after check_ipo_supported passes in Release builds.
# Leaving this ON would incorrectly set INTERPROCEDURAL_OPTIMIZATION on Debug targets.
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
        $<$<BOOL:${APPLE}>:-Wno-nullability-completeness>
        -pipe
        # -march=native / -mcpu=native tune for the machine being compiled on.
        # Always correct here since we build natively on each target machine.
        $<$<CONFIG:Release>:$<IF:$<BOOL:${APPLE}>,-mcpu=native,-march=native>>
        $<$<CONFIG:Release>:-O3>)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
endforeach()

# Explicit LTO flags per compiler (supplements INTERPROCEDURAL_OPTIMIZATION)
if(GCTA_IPO_ENABLED AND CMAKE_BUILD_TYPE MATCHES "^[Rr]elease$")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # AOCC only supports full LTO; upstream Clang/clang++ supports ThinLTO
        # which gives most of the benefit at a fraction of the link time.
        if(GCTA_COMPILER_IS_AOCC)
            target_compile_options(gcta64 PRIVATE -flto=full)
            target_link_options(gcta64    PRIVATE -flto=full)
        else()
            target_compile_options(gcta64 PRIVATE -flto=thin)
            target_link_options(gcta64    PRIVATE -flto=thin)
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # -flto=auto uses the same job count as make, no need to query core count
        target_compile_options(gcta64 PRIVATE -flto=auto)
        target_link_options(gcta64 PRIVATE -flto=auto)
    endif()
endif()

# Force DT_RPATH instead of DT_RUNPATH for all AOCC builds on Linux.
# AOCC/lld defaults to --enable-new-dtags which writes DT_RUNPATH.  Unlike
# DT_RPATH, DT_RUNPATH is searched *after* LD_LIBRARY_PATH — so if the HPC
# module system has prepended system lib dirs to LD_LIBRARY_PATH (which is
# very common), system libraries win even though the RPATH is correct.
# This applies to all AOCC-embedded RPATH entries: libomp, libstdc++, BLAS, etc.
# --disable-new-dtags reverts to DT_RPATH which is searched *before* LD_LIBRARY_PATH.
if(GCTA_COMPILER_IS_AOCC AND LINUX)
    target_link_options(gcta64 PRIVATE "LINKER:--disable-new-dtags")
    message(STATUS "gcta64: using DT_RPATH (--disable-new-dtags) for HPC library pinning")
endif()

# Static C++ runtime: embed libstdc++ and libgcc into the binary so it does not
# depend on the toolchain's shared .so files at runtime.  Useful on HPC clusters
# where the GCC module used at build time may not be available on compute nodes.
if(GCTA_STATIC_STDCXX AND LINUX)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_options(gcta64 PRIVATE -static-libstdc++ -static-libgcc)
        message(STATUS "gcta64: linking libstdc++ and libgcc statically")
    else()
        message(WARNING "GCTA_STATIC_STDCXX=ON is only supported with GCC or Clang/AOCC on Linux")
    endif()
endif()

add_library(Pgenlib
    third_party/Pgenlib/PgenReader.cpp
    third_party/plink-ng/2.0/include/pgenlib_ffi_support.cc
    third_party/plink-ng/2.0/include/pgenlib_misc.cc
    third_party/plink-ng/2.0/include/pgenlib_read.cc
    third_party/plink-ng/2.0/include/pgenlib_write.cc
    third_party/plink-ng/2.0/include/plink2_base.cc
    third_party/plink-ng/2.0/include/plink2_simd.cc
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
    # plink2_base.h requires FMA3 whenever 32-byte AVX2 float vectors are enabled.
    target_compile_options(Pgenlib PRIVATE -mavx2 -mfma -mbmi -mbmi2 -mlzcnt)
endif()

target_include_directories(Pgenlib PRIVATE ${GCTA_ZLIB_INCLUDE_DIRS})

if(BGEN_SUPPORT)
    target_link_libraries(gcta64 PRIVATE sqlite3)
    target_compile_definitions(gcta64 PRIVATE BGEN_SUPPORT)
endif()
target_link_libraries(gcta64 PRIVATE Pgenlib Threads::Threads)

# Link BLAS early so OpenBLAS/Accelerate/MKL cblas symbols take precedence
# over other cblas implementations in dyld (macOS) and ld.so (Linux) symbol resolution order
target_link_libraries(gcta64 PRIVATE ${GCTA_BLAS_LIBS})

# On Linux, re-link BLAS to ensure it appears last in NEEDED and
# its cblas symbols take precedence via ld.so's first-match ordering.
if(LINUX)
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

# Remove any duplicate entries added by compiler probes, BLAS detection, and
# toolchain logic (e.g. BLAS lib dir may match OpenMP dir on some installs).
list(REMOVE_DUPLICATES CMAKE_BUILD_RPATH)
list(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)

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
if(GCTA_COMPILER_IS_AOCC)
    set(GCTA_COMPILER_IS_AOCC_STR " (AOCC)")
else()
    set(GCTA_COMPILER_IS_AOCC_STR "")
endif()
message(STATUS "")
message(STATUS "==========================================")
message(STATUS "GCTAg Build Configuration Summary")
message(STATUS "==========================================")
message(STATUS "  Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "  BLAS/LAPACK:       ${GCTA_BLAS_PROVIDER}")
message(STATUS "  Compiler:          ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}${GCTA_COMPILER_IS_AOCC_STR} (${CMAKE_CXX_COMPILER})")
message(STATUS "  Archiver (AR):     ${CMAKE_AR}")
if(GCTA_COMPILER_IS_AOCC)
    message(STATUS "  GCC toolchain:     ${GCTA_GCC_TOOLCHAIN} (--gcc-toolchain)")
endif()
message(STATUS "  IPO/LTO:           ${GCTA_IPO_ENABLED}")
message(STATUS "  Install prefix:    ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  OpenMP:            ${OpenMP_CXX_VERSION}")
message(STATUS "  BGEN support:      ${BGEN_SUPPORT}")

message(STATUS "==========================================")
message(STATUS "")

set_property(DIRECTORY PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src"
    "${CMAKE_SOURCE_DIR}/main"
    "${CMAKE_SOURCE_DIR}/third_party/Pgenlib")