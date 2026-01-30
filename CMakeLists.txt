CMAKE_MINIMUM_REQUIRED(VERSION 3.18)

project(GCTA2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

# Find dependencies first before fetching libraries
find_package(Threads REQUIRED)

# GSL - Optional on all platforms
if(NOT WIN32)
    find_package(GSL)
    if(GSL_FOUND)
        message(STATUS "Found GSL: ${GSL_VERSION}")
    else()
        message(STATUS "GSL not found - some features may be disabled")
    endif()
endif()

#TODO: is this too complex?
# OpenMP - Cross-platform detection
find_package(OpenMP)
if(OpenMP_FOUND)
    message(STATUS "Found OpenMP: ${OpenMP_CXX_VERSION}")
else()
    if(APPLE)
        # Try Homebrew libomp on macOS
        find_path(OpenMP_C_INCLUDE_DIR omp.h
            HINTS /opt/homebrew/opt/libomp/include
                  /usr/local/opt/libomp/include)
        find_library(OpenMP_omp_LIBRARY omp
            HINTS /opt/homebrew/opt/libomp/lib
                  /usr/local/opt/libomp/lib)
        
        if(OpenMP_C_INCLUDE_DIR AND OpenMP_omp_LIBRARY)
            set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
            set(OpenMP_FOUND TRUE)
            
            add_library(OpenMP::OpenMP_C INTERFACE IMPORTED)
            set_target_properties(OpenMP::OpenMP_C PROPERTIES
                INTERFACE_COMPILE_OPTIONS "${OpenMP_C_FLAGS}"
                INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_C_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${OpenMP_omp_LIBRARY}")
            
            add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
            set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
                INTERFACE_COMPILE_OPTIONS "${OpenMP_CXX_FLAGS}"
                INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_C_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${OpenMP_omp_LIBRARY}")
            
            message(STATUS "Found Homebrew OpenMP: ${OpenMP_omp_LIBRARY}")
        endif()
    endif()
    
    if(NOT OpenMP_FOUND)
        message(FATAL_ERROR "OpenMP not found. Please install OpenMP.")
    endif()
endif()

# BLAS/LAPACK - Using CMake's built-in FindBLAS/FindLAPACK modules
option(USE_MKL "Prefer Intel MKL if available" ON)
set(USING_MKL FALSE)

# Set BLA_VENDOR to prefer specific implementations
if(USE_MKL AND (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i686"))
    set(BLA_VENDOR Intel10_64lp)  # Prefer MKL
elseif(APPLE)
    set(BLA_VENDOR OpenBLAS)  # Prefer OpenBLAS on macOS
else()
    set(BLA_VENDOR All)  # Accept any BLAS
endif()

# Add search paths for common installations
if(USE_MKL)
    list(APPEND CMAKE_PREFIX_PATH 
        "/opt/intel/mkl"
        "/opt/intel/oneapi/mkl/latest"
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
endif()

if(APPLE)
    # Add common OpenBLAS paths (Homebrew, conda, etc.)
    list(APPEND CMAKE_PREFIX_PATH 
        "/opt/homebrew/opt/openblas"
        "/usr/local/opt/openblas"
        "$ENV{CONDA_PREFIX}"  # Active conda environment
        "/Users/alexleonard/Documents/Tiergenomik/software/miniconda3/envs/rust")
endif()

find_package(BLAS)
find_package(LAPACK)

if(BLAS_FOUND)
    message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
    # Detect if we're using MKL by checking library names
    string(FIND "${BLAS_LIBRARIES}" "mkl" MKL_POS)
    if(MKL_POS GREATER -1)
        set(USING_MKL TRUE)
        message(STATUS "Using Intel MKL")
    endif()
    
    # Detect if we're using OpenBLAS
    string(FIND "${BLAS_LIBRARIES}" "openblas" OPENBLAS_POS)
    if(OPENBLAS_POS GREATER -1)
        set(USING_OPENBLAS TRUE)
        message(STATUS "Using OpenBLAS")
    endif()
else()
    message(WARNING "No BLAS found. Performance may be degraded.")
endif()

if(LAPACK_FOUND)
    message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
else()
    message(WARNING "No LAPACK found. Performance may be degraded.")
endif()

# Get include directories for vecLib on macOS (only if using Accelerate, not OpenBLAS)
if(APPLE AND BLAS_FOUND AND NOT USING_OPENBLAS)
    execute_process(
        COMMAND xcrun --sdk macosx --show-sdk-path
        OUTPUT_VARIABLE MACOSX_SDK_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(VECLIB_INCLUDE_DIR "${MACOSX_SDK_PATH}/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers")
    set(USING_ACCELERATE TRUE)
    message(STATUS "Using Accelerate framework")
endif()

# Fetch 3rd parties using FetchContent
include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# zlib
FetchContent_Declare(
    zlib
    URL https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(ZLIB_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(zlib)

# Eigen - Updated to 5.0 (header-only)
# Disable Pardiso globally if not using MKL (MUST be set before fetching Eigen)
if(NOT USING_MKL)
    add_compile_definitions(EIGEN_NO_PARDISO)
    message(STATUS "Setting EIGEN_NO_PARDISO globally (no MKL detected)")
endif()

FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_TEST_MKL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(Eigen3)

# Spectra - Updated to 1.2 (header-only)
FetchContent_Declare(
    Spectra
    URL https://github.com/yixuan/spectra/archive/refs/tags/v1.2.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(CMAKE_SKIP_INSTALL_RULES ON)
FetchContent_MakeAvailable(Spectra)
set(CMAKE_SKIP_INSTALL_RULES OFF)

# zstd
FetchContent_Declare(
    zstd
    URL https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR build/cmake
)
set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
set(ZSTD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ZSTD_BUILD_SHARED ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(zstd)

FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.90.0/boost-1.90.0-cmake.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(BOOST_ENABLE_CMAKE ON CACHE BOOL "" FORCE)
set(BOOST_INCLUDE_LIBRARIES "algorithm;math;crc" CACHE STRING "" FORCE)
FetchContent_MakeAvailable(Boost)

# SQLite
FetchContent_Declare(
    sqlite
    URL https://www.sqlite.org/2025/sqlite-amalgamation-3480000.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(sqlite)
if(NOT TARGET sqlite3)
    add_library(sqlite3 STATIC ${sqlite_SOURCE_DIR}/sqlite3.c)
    target_include_directories(sqlite3 PUBLIC ${sqlite_SOURCE_DIR})
    target_compile_definitions(sqlite3 PRIVATE 
        SQLITE_ENABLE_RTREE=1
        SQLITE_ENABLE_FTS5=1)
endif()

# Create interface library for common dependencies
add_library(gcta_common INTERFACE)
target_include_directories(gcta_common INTERFACE
    ${CMAKE_SOURCE_DIR}/include
    "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/include"
    "${CMAKE_SOURCE_DIR}/third_party/plink-ng/2.0/simde"
    ${zstd_SOURCE_DIR}/lib
    ${zlib_SOURCE_DIR}
    ${zlib_BINARY_DIR})

if(USING_MKL AND MKL_INCLUDE)
    target_include_directories(gcta_common INTERFACE ${MKL_INCLUDE})
    # Configure Eigen with full MKL support
    target_compile_definitions(gcta_common INTERFACE 
        EIGEN_USE_MKL_ALL
        EIGEN_USE_BLAS
        EIGEN_USE_LAPACKE
        EIGEN_USE_PARDISO)
    message(STATUS "Eigen configured with MKL support")
elseif(USING_OPENBLAS)
    # OpenBLAS: find and add include directories for cblas.h and lapacke.h
    get_filename_component(OPENBLAS_LIB_DIR "${BLAS_LIBRARIES}" DIRECTORY)
    get_filename_component(OPENBLAS_PREFIX "${OPENBLAS_LIB_DIR}" DIRECTORY)
    find_path(OPENBLAS_INCLUDE_DIR cblas.h
        HINTS "${OPENBLAS_PREFIX}/include"
              "${OPENBLAS_PREFIX}/include/openblas"
              /opt/homebrew/opt/openblas/include
              /usr/local/opt/openblas/include
        NO_DEFAULT_PATH)
    
    # If not found, try default paths
    if(NOT OPENBLAS_INCLUDE_DIR)
        find_path(OPENBLAS_INCLUDE_DIR cblas.h
            HINTS "${OPENBLAS_PREFIX}/include"
                  "${OPENBLAS_PREFIX}/include/openblas"
                  /opt/homebrew/opt/openblas/include
                  /usr/local/opt/openblas/include)
    endif()
    
    if(OPENBLAS_INCLUDE_DIR)
        # Use BEFORE to prioritize OpenBLAS headers over system headers
        target_include_directories(gcta_common BEFORE INTERFACE ${OPENBLAS_INCLUDE_DIR})
        message(STATUS "Found OpenBLAS headers: ${OPENBLAS_INCLUDE_DIR}")
    endif()
    
    # Configure Eigen and plink-ng to use OpenBLAS
    target_compile_definitions(gcta_common INTERFACE 
        EIGEN_USE_BLAS
        EIGEN_USE_LAPACKE
        USE_OPENBLAS  # Tell plink-ng to use OpenBLAS instead of Accelerate
        __ACCELERATE__=0)  # Prevent Accelerate framework inclusion
    
    # Add compiler flag to exclude Accelerate framework on macOS
    if(APPLE)
        target_compile_options(gcta_common INTERFACE 
            -DEIGEN_NO_AUTOMATIC_LINKING
            -Wno-deprecated-declarations)
    endif()
    
    message(STATUS "Eigen configured with OpenBLAS + LAPACKE")
elseif(USING_ACCELERATE)
    target_include_directories(gcta_common INTERFACE ${VECLIB_INCLUDE_DIR})
    # Configure Eigen with Accelerate framework (EIGEN_NO_PARDISO set globally above)
    target_compile_definitions(gcta_common INTERFACE 
        EIGEN_USE_BLAS
        EIGEN_USE_LAPACKE)
    message(STATUS "Eigen configured with Accelerate framework")
else()
    # Generic BLAS/LAPACK (EIGEN_NO_PARDISO set globally above)
    target_compile_definitions(gcta_common INTERFACE EIGEN_USE_BLAS)
    message(STATUS "Eigen configured with generic BLAS/LAPACK")
endif()

if(GSL_FOUND)
    target_include_directories(gcta_common INTERFACE ${GSL_INCLUDE_DIRS})
endif()

target_link_libraries(gcta_common INTERFACE
    Boost::algorithm Boost::math Boost::crc
    Eigen3::Eigen
    Spectra
    sqlite3)

if(OpenMP_CXX_FOUND)
    target_link_libraries(gcta_common INTERFACE OpenMP::OpenMP_CXX)
endif()

if(BLAS_FOUND AND LAPACK_FOUND)
    target_link_libraries(gcta_common INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# Build individual library targets from source files
file(GLOB SRCS "${CMAKE_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRCS "${CMAKE_SOURCE_DIR}/src/main.cpp")
list(APPEND SRCS "${CMAKE_SOURCE_DIR}/src/acat.c")

set(GCTA_LIBS "")
foreach(lib_src ${SRCS})
    get_filename_component(lib_name "${lib_src}" NAME_WE)
    string(TOLOWER ${lib_name} lib_name)
    add_library(${lib_name} "${lib_src}")
    target_link_libraries(${lib_name} PUBLIC gcta_common)
    list(APPEND GCTA_LIBS ${lib_name})
endforeach()

# Build main executable
add_executable(gcta64 "${CMAKE_SOURCE_DIR}/src/main.cpp")

add_subdirectory(main)
target_link_libraries(mainV1 PUBLIC gcta_common)

add_subdirectory(third_party)
target_include_directories(Pgenlib PRIVATE ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})

# Link main executable
target_link_libraries(gcta64
    gcta_common
    mainV1
    ${GCTA_LIBS}
    Pgenlib
    libzstd
    Threads::Threads
    zlib)

if(GSL_FOUND)
    target_link_libraries(gcta64 GSL::gsl GSL::gslcblas)
endif()

# Copy compile_commands.json to source directory for IDE support
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
    file(COPY "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json"
         DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# Installation
if(LINUX OR APPLE)
    set(CMAKE_INSTALL_RPATH "$<IF:$<PLATFORM_ID:Darwin>,@executable_path/../lib,$ORIGIN/../lib>")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)
    
    install(TARGETS gcta64 mainV1 ${GCTA_LIBS}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
    
    if(LINUX)
        # Linux desktop integration files
        install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.desktop
            DESTINATION share/applications)
        install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.png
            DESTINATION share/icons/hicolor/256x256/apps)
        install(FILES ${CMAKE_SOURCE_DIR}/resources/gcta64.svg
            DESTINATION share/icons/hicolor/scalable/apps)
    endif()
elseif(WIN32)
    install(TARGETS gcta64 mainV1 ${GCTA_LIBS} DESTINATION bin)

    file(GLOB RUNTIME_DLLS "${CMAKE_BINARY_DIR}/*.dll")
    set(INSTALL_DLLS ${RUNTIME_DLLS} ${zlib_BINARY_DIR}/zlib.dll)
    
    # Intel MKL DLLs if using MKL
    if(USING_MKL AND MKL_ROOT)
        list(APPEND INSTALL_DLLS
            ${MKL_ROOT}/bin/mkl_intel_thread.2.dll
            ${MKL_ROOT}/bin/mkl_avx2.2.dll
            ${MKL_ROOT}/bin/mkl_core.2.dll
            ${MKL_ROOT}/../../compiler/latest/bin/libiomp5md.dll)
    endif()
    
    install(FILES ${INSTALL_DLLS} DESTINATION bin)
endif()
