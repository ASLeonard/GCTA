# GCTAg

`GCTAg` is a fork of `GCTA`, aiming at both scaling to larger datasets with improved compute, but also exploiting the typically dense GRMs of agricultural (hence GCT**Ag**) datasets.

A summary of the new flags is given in [here](docs/changes/GCTAg.md), with some (incomplete) notes on the linear algebra optimisations [here](docs/development/linalg_optimizations.md).
A manuscript summarising this work is in preperation.

## GCTA

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

Alex Leonard programmed the `GCTAg` fork, improving performance in existing `GCTA` bottlenecks as well as adding additional features (Hutch++, Woodbury, etc) to scale to larger cohorts.

Contributions to the development of the methods implemented in GCTA (e.g., GREML methods, COJO, mtCOJO, MLMA-LOCO, fastBAT, fastGWA and fastGWA-GLMM) can be found in the corresponding publications (https://yanglab.westlake.edu.cn/software/gcta/index.html#Overview).


## Questions and Help Requests
If you have any bug reports or questions please send an email to Jian Yang at <jian.yang@westlake.edu.cn>.


## Compilation

#### Requirements

x86\_64 and ARM CPUs are supported.
2. A BLAS backend: [AOCL](https://www.amd.com/en/developer/aocl.html) (recommended on AMD EPYC/Zen4), [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) 2017 or above, or OpenBLAS. The backend is no longer auto-detected; supply it explicitly via the `GCTA_BLAS_LIBRARY` and `GCTA_BLAS_INCLUDE_DIR` CMake variables (see Build below).
3. Eigen 5.0.1 (downloaded automatically via CMake FetchContent; the old 3.3.7 pin is no longer required).
4. CMake >= 3.28
5. BOOST >= 1.90 (Boost.Math is used for statistical distributions; downloaded via FetchContent, with its zlib dependency resolved automatically)
6. zlib >= 1.3.2
8. zstd >= 1.5.7
9. [Spectra](https://spectralib.org/) >= 1.2.0

Optional for BGEN support
1. sqlite3 >= 3.31.1

#### Linux

1. Kernel version >= 2.6.28.
2. GCC version >= 13 or Clang >= 17, with C++23 support. An AOCC toolchain file is available under `cmake/` for building with AMD's compiler.

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

Point CMake at your BLAS backend explicitly, e.g. for AOCL on the ETH Euler/Leonhard clusters:

```sh
cmake -DCMAKE_BUILD_TYPE=Release \
      -DGCTA_BLAS_LIBRARY=/path/to/aocl/lib/libblis-mt.so \
      -DGCTA_BLAS_INCLUDE_DIR=/path/to/aocl/include \
      -G Ninja -B build/Release -S .
```

RPATH is prepended automatically so the built binary resolves the configured BLAS library at runtime rather than falling back to a system default.

On Windows, you should use the toolchain file in `cmake/win-toolchain.cmake`:

``` sh
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="cmake/win-toolchain.cmake" -G Ninja -B build/Release -S .
```

##### Compile

```sh
cmake --build build/Release
```

The executable binary will be generated under `build/Release`.