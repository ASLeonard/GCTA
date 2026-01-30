# CMake Build System Improvements

## Overview
This document describes the modernization and improvements made to the GCTA CMake build system.

## Key Improvements

### 1. **Flexible BLAS/LAPACK Detection**
The build system now intelligently detects and uses the best available linear algebra library:

- **Intel MKL** (preferred on x86_64 platforms):
  - Automatically detected via CMake's `find_package(MKL CONFIG)`
  - Searches common installation paths on Linux and Windows
  - Falls back to alternatives if not found
  - Can be disabled with `-DUSE_MKL=OFF`

- **Apple Accelerate** (macOS):
  - Uses the built-in Accelerate framework
  - Includes vecLib headers for BLAS/LAPACK
  - No additional installation required

- **Generic BLAS/LAPACK** (fallback):
  - Works with OpenBLAS, ATLAS, or system BLAS/LAPACK
  - Automatically detected on Linux systems

### 2. **Updated Dependencies**

| Library | Old Version | New Version |
|---------|-------------|-------------|
| Boost   | 1.86.0      | 1.90.0      |
| Eigen   | 3.4.0       | 5.0.0       |
| Spectra | 1.1.0       | 1.2.0       |
| SQLite  | 3.46.0      | 3.47.2      |
| zlib    | 1.3.1       | 1.3.1       |
| zstd    | 1.5.6       | 1.5.6       |

### 3. **Improved OpenMP Detection**

- Cross-platform OpenMP detection
- Automatic fallback to Homebrew libomp on macOS
- Proper error messages if OpenMP is not found
- Creates imported targets for consistent linking

### 4. **Better GSL Handling**

- GSL is now optional on all platforms (not just non-Windows)
- Graceful degradation if GSL is not available
- Clear status messages about GSL availability

### 5. **Modern CMake Practices**

- Minimum CMake version raised to 3.18
- Position-independent code enabled by default
- Standard requirements enforced with `CMAKE_CXX_STANDARD_REQUIRED`
- Proper use of imported targets throughout
- Better dependency propagation

## Building the Project

### Prerequisites

#### Common Requirements (All Platforms)
- CMake 3.18 or later
- C++17 compatible compiler (GCC 6.1+, Clang 8.0+)
- OpenMP support

#### Platform-Specific Requirements

##### macOS
```bash
# Install required dependencies via Homebrew
brew install cmake libomp

# Optional but recommended
brew install gsl openblas
```

##### Linux
```bash
# Ubuntu/Debian
sudo apt-get install cmake build-essential libgomp1

# Optional: Install MKL for better performance
# Download from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

# Or use OpenBLAS as alternative
sudo apt-get install libopenblas-dev liblapack-dev

# Optional: GSL
sudo apt-get install libgsl-dev
```

##### Windows
- Visual Studio 2019 or later (with C++ support)
- LLVM/Clang for Windows
- OpenMP runtime
- Optional: Intel MKL

### Build Instructions

#### Standard Build
```bash
cd GCTA
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

#### Build with Specific BLAS/LAPACK

**Force MKL usage:**
```bash
cmake -DUSE_MKL=ON -DMKL_DIR=/path/to/mkl ..
```

**Disable MKL (use system BLAS/LAPACK):**
```bash
cmake -DUSE_MKL=OFF ..
```

#### Build Types
```bash
# Release build (default, optimized)
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug build (with symbols)
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Release with debug info
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
```

#### Custom Compiler
```bash
# Use specific compiler
cmake -DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11 ..

# Or with Clang
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
```

### Installation
```bash
# Install to system (requires sudo on Linux/macOS)
sudo cmake --install . --prefix /usr/local

# Install to custom location
cmake --install . --prefix ~/myapps/gcta
```

## Configuration Options

| Option | Default | Description |
|--------|---------|-------------|
| `USE_MKL` | ON | Use Intel MKL if available |
| `CMAKE_BUILD_TYPE` | Release | Build type (Release/Debug/RelWithDebInfo) |
| `MKL_DIR` | (auto) | Path to MKL installation |

## Troubleshooting

### MKL Not Found
If you have MKL installed but CMake can't find it:
```bash
cmake -DMKL_DIR=/opt/intel/oneapi/mkl/latest/lib/cmake/mkl ..
```

Or disable MKL and use alternatives:
```bash
cmake -DUSE_MKL=OFF ..
```

### OpenMP Not Found (macOS)
```bash
# Install via Homebrew
brew install libomp

# If still not found, specify manually
cmake -DOpenMP_C_INCLUDE_DIR=/opt/homebrew/opt/libomp/include \
      -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib ..
```

### Missing Dependencies
The build system automatically fetches most dependencies via CPM. If network issues occur:
```bash
# Set CPM cache to reuse downloads
export CPM_SOURCE_CACHE=$HOME/.cache/CPM
cmake ..
```

## Performance Considerations

### BLAS/LAPACK Performance Hierarchy
1. **Intel MKL** - Fastest on Intel CPUs (x86_64)
2. **Apple Accelerate** - Optimized for Apple Silicon and Intel Macs
3. **OpenBLAS** - Good general-purpose alternative
4. **Generic BLAS/LAPACK** - Slower but portable

### Recommended Configurations

**For maximum performance on Linux (x86_64):**
```bash
# Install Intel MKL
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

cmake -DUSE_MKL=ON -DCMAKE_BUILD_TYPE=Release ..
```

**For Apple Silicon:**
```bash
# Uses Accelerate framework automatically
cmake -DCMAKE_BUILD_TYPE=Release ..
```

**For generic Linux without MKL:**
```bash
# Install OpenBLAS
sudo apt-get install libopenblas-dev

cmake -DUSE_MKL=OFF -DCMAKE_BUILD_TYPE=Release ..
```

## Migration from Old Build System

### Key Changes
1. MKL is no longer required - builds work with any BLAS/LAPACK
2. Hardcoded paths removed - everything auto-detected
3. Library versions updated to latest stable releases
4. Better cross-platform support

### Old Paths No Longer Needed
These environment variables/paths are no longer required:
- `MKL_DIR` (auto-detected, but can still be set)
- `OMP_LIBRARY` (auto-detected)
- Hardcoded `/opt/intel/mkl` paths
- Hardcoded Homebrew paths (still checked as fallbacks)

## Testing the Build

After building, verify the executable:
```bash
./gcta64 --version
./gcta64 --help
```

## Contributing

When modifying the build system:
- Test on multiple platforms (Linux, macOS, Windows)
- Test with and without MKL
- Test with different compilers (GCC, Clang)
- Verify dependency updates don't break functionality
- Update this documentation

## Support

For build issues:
1. Check this documentation first
2. Verify all prerequisites are installed
3. Try the troubleshooting steps above
4. Check CMake output for specific error messages
5. File an issue with complete build log
