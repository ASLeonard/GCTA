CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

project(GCTA2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_C_STANDARD 23)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (CMAKE_COMPILER_IS_GNUCC)
    MESSAGE("Compiling by GCC ${CMAKE_CXX_COMPILER_VERSION}")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 6.0))
        MESSAGE(FATAL_ERROR "GCTA can't build without GCC version >= 6.1")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    MESSAGE("Compiling by Clang ${CMAKE_CXX_COMPILER_VERSION}")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0))
        MESSAGE(FATAL_ERROR "GCTA can't build without Clang version >= 8.0")
    endif()

    # If using Homebrew LLVM on macOS, ensure we use its libc++
    if(APPLE AND CMAKE_CXX_COMPILER MATCHES "/opt/homebrew/opt/llvm")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        link_directories(/opt/homebrew/opt/llvm/lib/c++)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/opt/homebrew/opt/llvm/lib/c++")
        list(APPEND CMAKE_BUILD_RPATH "/opt/homebrew/opt/llvm/lib/c++")
        message(STATUS "Using Homebrew LLVM libc++ from /opt/homebrew/opt/llvm/lib/c++")
    endif()
else()
    MESSAGE(FATAL_ERROR "GCTA can only compile with GCC or Clang")
endif()

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
    message("Setting to release build")
endif()

# Set optimization flags for Release builds
if(CMAKE_BUILD_TYPE MATCHES "^[Rr]elease$")
    if(APPLE)
        # Apple Silicon optimizations
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mcpu=native -pipe ${CMAKE_CXX_FLAGS_RELEASE}")
        set(CMAKE_C_FLAGS_RELEASE "-O3 -mcpu=native -pipe ${CMAKE_C_FLAGS_RELEASE}")
    else()
        # x86_64 optimizations
        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -pipe ${CMAKE_CXX_FLAGS_RELEASE}")
        set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -pipe ${CMAKE_C_FLAGS_RELEASE}")
    endif()
    message(STATUS "Release optimization flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

option(GCTA_USE_MKL "Use Intel MKL for BLAS/LAPACK (non-Apple only)" ON)
option(GCTA_LIST_TARGETS "List imported and local targets during configure" OFF)
option(GCTA_UNITY_BUILD "Enable unity builds for faster compilation" OFF)
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
set(CMAKE_DISABLE_FIND_PACKAGE_CLANG_FORMAT ON CACHE BOOL "" FORCE)

set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
set(ZSTD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ZSTD_BUILD_SHARED ON CACHE BOOL "" FORCE)
set(BOOST_ENABLE_CMAKE ON CACHE BOOL "" FORCE)
set(BOOST_INCLUDE_LIBRARIES "algorithm;math;crc" CACHE STRING "" FORCE)

FetchContent_Declare(
    zlib
    URL https://github.com/madler/zlib/archive/refs/tags/v1.3.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
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
    SOURCE_SUBDIR "build/cmake")
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-cmake.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(
    sqlite
    URL https://www.sqlite.org/2026/sqlite-amalgamation-3510200.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)

set(SPECTRA_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(CMAKE_SKIP_INSTALL_RULES ON)

FetchContent_MakeAvailable(zlib eigen spectra zstd Boost sqlite)

set(CMAKE_SKIP_INSTALL_RULES OFF)
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

set(OPENBLAS_INCLUDE_DIR "" CACHE PATH "Path to OpenBLAS headers when GCTA_USE_MKL=OFF")
set(GCTA_BLAS_LIBS "")
set(GCTA_BLAS_INCLUDE_DIRS "")
set(GCTA_BLAS_DEFINES "")

if(APPLE)
    if(NOT CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        message(FATAL_ERROR "GCTA macOS support is arm64 only.")
    endif()
    set(GCTA_USE_MKL OFF)
    find_library(ACCELERATE_LIBRARY Accelerate REQUIRED)
    list(APPEND GCTA_BLAS_LIBS ${ACCELERATE_LIBRARY})
    list(APPEND GCTA_BLAS_DEFINES GCTA_USE_ACCELERATE ACCELERATE_NEW_LAPACK ACCELERATE_LAPACK_ILP64)
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
else()
    if(GCTA_USE_MKL)
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
        list(APPEND GCTA_BLAS_DEFINES GCTA_USE_MKL)
        set(GCTA_BLAS_PROVIDER "MKL" CACHE INTERNAL "")
    else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        list(APPEND GCTA_BLAS_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        list(APPEND GCTA_BLAS_DEFINES GCTA_USE_OPENBLAS)
        set(GCTA_BLAS_PROVIDER "BLAS/LAPACK" CACHE INTERNAL "")
        if(OPENBLAS_INCLUDE_DIR)
            list(APPEND GCTA_BLAS_INCLUDE_DIRS ${OPENBLAS_INCLUDE_DIR})
        endif()
    endif()
endif()

message(STATUS "BLAS/LAPACK provider: ${GCTA_BLAS_PROVIDER}")

set(COMMON_INCLUDES ${CMAKE_SOURCE_DIR}/include)
set(MAIN_SOURCE "${CMAKE_SOURCE_DIR}/src/main.cpp")

file(GLOB SRCS "${CMAKE_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${MAIN_SOURCE}")
list(APPEND SRCS "${CMAKE_SOURCE_DIR}/src/acat.c")

set(GCTA_LIBS "")
foreach(lib_src ${SRCS})
    get_filename_component(file_name "${lib_src}" NAME_WE)
    string(TOLOWER ${file_name} lib_name)
    add_library(${lib_name} "${lib_src}")

    target_include_directories(${lib_name} PUBLIC
        ${COMMON_INCLUDES}
        ${GCTA_BLAS_INCLUDE_DIRS}
        # PLINK2.0 pgen library;
        "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include"
        "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/simde")

    target_link_libraries(${lib_name}
        Boost::algorithm Boost::math Boost::crc
        Eigen3::Eigen
        Spectra
        sqlite3)

    list(APPEND GCTA_LIBS ${lib_name})
endforeach(lib_src)

add_executable(gcta64 ${MAIN_SOURCE})

add_subdirectory(main)
target_include_directories(mainV1 PUBLIC ${COMMON_INCLUDES} ${GCTA_BLAS_INCLUDE_DIRS})
target_link_libraries(mainV1 Eigen3::Eigen)

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
endforeach()

# Apply common includes and definitions to all targets
foreach(target IN LISTS GCTA_TARGETS)
    target_include_directories(${target} PRIVATE
        ${zstd_SOURCE_DIR}/lib
        ${zlib_SOURCE_DIR}
        ${zlib_BINARY_DIR}
        ${GCTA_BLAS_INCLUDE_DIRS})
    target_compile_definitions(${target} PRIVATE ${GCTA_BLAS_DEFINES})
    if(APPLE)
        target_compile_options(${target} PRIVATE -Wno-nullability-completeness)
    endif()
    if(NOT WIN32)
        target_include_directories(${target} PRIVATE ${GSL_INCLUDE_DIRS})
    endif()
endforeach()

# Link OpenMP to all targets
foreach(target IN LISTS GCTA_TARGETS)
    target_link_libraries(${target} OpenMP::OpenMP_CXX)
endforeach()

add_subdirectory(third_party)

target_include_directories(gcta64 PRIVATE
    ${COMMON_INCLUDES}
    ${GCTA_BLAS_INCLUDE_DIRS})

target_include_directories(Pgenlib PRIVATE ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})

target_link_libraries(gcta64 mainV1 ${GCTA_LIBS} Pgenlib sqlite3 libzstd Threads::Threads zlibstatic)

if (NOT WIN32)
    target_link_libraries(gcta64 GSL::gsl GSL::gslcblas)
endif()

foreach(target IN LISTS GCTA_TARGETS)
    target_link_libraries(${target} ${GCTA_BLAS_LIBS})
endforeach()

# Testing has some problems currently
#enable_testing()
#ADD_SUBDIRECTORY(test)

IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
    )
ENDIF()

# List all targets (verbose, only if requested)
if(CMAKE_VERBOSE_MAKEFILE OR GCTA_LIST_TARGETS)
    get_property(imported_targets DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
    foreach(target IN LISTS imported_targets)
        message(STATUS "imported target: ${target}")
    endforeach()

    get_property(local_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
    foreach(t IN LISTS local_targets)
        message(STATUS "built target: ${t}")
    endforeach()
endif()

# Installation
if(LINUX)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

    # Targets
    install(TARGETS gcta64 mainV1 ${GCTA_LIBS}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

    # Resource
    install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.desktop
        DESTINATION share/applications)
    install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.png
        DESTINATION share/icons/hicolor/256x256/apps)
    install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.svg
        DESTINATION share/icons/hicolor/scalable/apps)
elseif(WIN32)
    install(TARGETS gcta64 mainV1 ${GCTA_LIBS} DESTINATION bin)

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
elseif(APPLE)
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

    # Install target executables and libraries
    install(TARGETS gcta64 mainV1 ${GCTA_LIBS}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

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
