# GCTA
GCTA (Genome-wide Complex Trait Analysis) is a software package, which was initially developed to estimate the proportion of phenotypic variance explained by all genome-wide SNPs for a complex trait but has been extensively extended for many other analyses of data from genome-wide association studies (GWASs). Please see the software website through the link below for more information.

Software website: https://yanglab.westlake.edu.cn/software/gcta/
License: GPLv3 (some parts of the code are released under LGPL as detailed in the files).


## Credits  
Jian Yang developed the original version (before v1.90) of the software (with supports from Peter Visscher, Mike Goddard and Hong Lee) and currently maintains the software.

Zhili Zheng programmed the fastGWA, fastGWA-GLMM and fastGWA-BB modules, rewrote the I/O and GRM modules, improved the GREML and bivariate GREML modules, extended the PCA module, and improved the SBLUP module.  

Zhihong Zhu programmed the mtCOJO and GSMR modules and improved the COJO module.  

Longda Jiang and Hailing Fang developed the ACAT-V module.  

Jian Zeng rewrote the GCTA-HEreg module.  

Andrew Bakshi contributed to the GCTA-fastBAT module.

Angli Xue improved the GSMR module.

Robert Maier improved the GCTA-SBLUP module.

Wujuan Zhong and Judong Shen programmed the fastGWA-GE module. 

Contributions to the development of the methods implemented in GCTA (e.g., GREML methods, COJO, mtCOJO, MLMA-LOCO, fastBAT, fastGWA and fastGWA-GLMM) can be found in the corresponding publications (https://yanglab.westlake.edu.cn/software/gcta/index.html#Overview).


## Questions and Help Requests
If you have any bug reports or questions please send an email to Jian Yang at <jian.yang@westlake.edu.cn>.


## Compilation

### Quick Start

**See [docs/build/QUICK_START.md](docs/build/QUICK_START.md) for platform-specific instructions.**

For detailed build information and troubleshooting, see [BUILD_IMPROVEMENTS.md](BUILD_IMPROVEMENTS.md).

#### Requirements

**Core Requirements (all platforms):**
- CMake >= 3.18
- C++17 compatible compiler (GCC 6.1+, Clang 8.0+, or MSVC)
- OpenMP support

**Optional (Auto-downloaded by CMake):**
- Boost 1.90.0 (algorithm, math, crc)
- Eigen 5.0.0
- Spectra 1.2.0
- zlib 1.3.1
- zstd 1.5.6
- sqlite3 3.47.2

**Linear Algebra Library (one of):**
- Intel MKL 2017+ (preferred on x86\_64, optional)
- Apple Accelerate (macOS, automatic)
- OpenBLAS (Linux/cross-platform alternative)
- Generic BLAS/LAPACK (universal fallback)

**Platform-Specific:**
- GNU Scientific Library (gsl) - optional but recommended

#### Key Improvements in Build System

✅ **MKL optional** - builds work without MKL using OpenBLAS or system BLAS/LAPACK  
✅ **Auto-detection** - no hardcoded paths, works out of the box on most systems  
✅ **Updated dependencies** - Boost 1.90, Eigen 5.0, Spectra 1.2  
✅ **Better error messages** - clear indication of what's missing  
✅ **Cross-platform** - improved support for ARM64, macOS, and Linux  

#### Simple Build

```bash
cd GCTA
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j$(nproc)
./gcta64 --version
```

#### Build Without MKL

If you don't have Intel MKL or prefer to use OpenBLAS:

```bash
cmake -DUSE_MKL=OFF -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j$(nproc)
```

#### Linux

1. Kernel version >= 2.6.28 (otherwise the Intel MKL library doesn't work).
2. GCC version >= 6.1 with C++ 11 support.

#### Before compilation 

Update [plink_ng](https://github.com/chrchang/plink-ng) submodule first.

```sh
git submodule update --init
```

On Windows, apply the patch under the `third_party` directory to the `plink-ng`.

#### Build

##### CMake Configuration

On MacOS and Linux, use following command to generate the build system:

```sh
cmake -DCMAKE_BUILD_TYPE=Release -G Ninja -B build/Release -S .
```

On Windows, you should use the toolchain file in `cmake/win-toolchain.cmake`:

``` sh
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="cmake/win-toolchain.cmake" -G Ninja -B build/Release -S .
```

##### Compile

```sh
cmake --build build/Release
```

The executable binary will be generated under `build/Release`.

