# Linear Algebra Optimisations in the GCTA Refactor

This document catalogues performance-relevant changes introduced across the GCTA codebase. It is split into two parts:

- **Part I — Core Performance** covers algorithmic and numerical improvements: new operations, loop restructuring, memory layout, and a stochastic trace estimator. Mathematical context is given for each change.
- **Part II — Supporting Changes** covers portability (MKL removal, platform guards), correctness (Boost.Math), I/O, and new CLI flags that unlock the Part I features.

---

# Part I: Core Performance Improvements

---

## 1. GRM Construction (`src/GRM.cpp`)

The genomic relationship matrix is defined as:

$$G = \frac{1}{m} Z Z^T, \qquad Z_{ij} = \frac{x_{ij} - 2p_j}{\sqrt{2p_j(1-p_j)}}$$

where $x_{ij}$ is the raw allele count for individual $i$ at marker $j$ and $p_j$ is the allele frequency. The computation proceeds in blocks of $k$ markers at a time, applying rank-$k$ updates $G \mathrel{+}= \frac{1}{m} Z_k Z_k^T$ per block.

### 1.1 Raw BLAS → Eigen `selfadjointView` / `rankUpdate`

**Before**
```cpp
// #if guards for x86 vs non-x86 ABI differences
#if GCTA_CPU_x86
    dsyrk(&uplo, &notrans, &n, &curNumValidMarkers, &alpha, stdGeno, &n_sample, &beta, grm, &m);
    dgemm(&notrans, &trans, &m, &s_n, &curNumValidMarkers, &alpha,
          stdGeno + part_keep_indices.first, &n_sample, stdGeno, &n_sample, &beta, grm, &m);
#else
    dsyrk_(/* Fortran ABI */);
    dgemm_(/* Fortran ABI */);
#endif
```

**After**
```cpp
// Full-sample partition (symmetric rank-k update = DSYRK)
Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>> A(stdGeno, grm_n, k, Eigen::OuterStride<>(stdGenoLD));
Eigen::Map<Eigen::MatrixXd>(grm, grm_n, grm_n)
    .selfadjointView<Eigen::Lower>().rankUpdate(A, alpha);

// Cross-partition (rectangular block = DGEMM + DSYRK)
Eigen::Map<...> A_top(stdGeno, grm_s_n, k, stride);
Eigen::Map<...> A_bot(stdGeno + part_keep_indices.first, grm_m, k, stride);
Eigen::Map<Eigen::MatrixXd>(grm, grm_m, grm_s_n).noalias() += alpha * (A_bot * A_top.transpose());
Eigen::Map<Eigen::MatrixXd>(grm_start, grm_m, grm_m)
    .selfadjointView<Eigen::Lower>().rankUpdate(A_bot, alpha);
```

`selfadjointView<Lower>::rankUpdate` maps to BLAS `DSYRK`, writing only the lower triangle. This removes platform-specific `#if` guards and makes the intent explicit to both the compiler and the reader.

### 1.2 Aligned column stride (`stdGenoLD`)

AVX/AVX-512 BLAS kernels require 64-byte aligned leading dimensions. The stride is now rounded up:

$$\text{stdGenoLD} = \left\lceil \frac{n}{8} \right\rceil \times 8$$

```cpp
stdGenoLD = (grm_n + 7) / 8 * 8;
```

Every column of the standardised genotype buffer is now 64-byte aligned, satisfying the alignment requirement of AVX/AVX-512 BLAS kernels and eliminating the scalar fallback path that fires for unaligned loads.

### 1.3 Caching derived constants

The old code used `static int` local variables inside `calculate_GRM_blas` that were re-initialised every call (harmless correctness-wise but confusing and thread-unsafe across calls from different `GRM` instances). These are now computed once in the constructor:

```cpp
grm_n = static_cast<int>(part_keep_indices.second) + 1;
grm_m = static_cast<int>(part_keep_indices.second) - static_cast<int>(part_keep_indices.first) + 1;
grm_s_n = grm_n - grm_m;
grm_bytes_std_geno = sizeof(double) * grm_n;
stdGenoLD = (grm_n + 7) / 8 * 8;
```

### 1.4 Loop hoisting: reducing OpenMP fork/join barriers

**Before** — two parallel regions *inside* the `numNblock` loop:
```
for i in 0..numNblock:
    #pragma omp parallel for   ← fork/join × numNblock
    fill sample_miss[i]
    #pragma omp parallel for   ← fork/join × numNblock
    N_thread(sample_miss[i])
delete[] sample_miss
```

**After** — fill all blocks first, then one parallel region over pairs:
```
for i in 0..numNblock:
    fill sampleMissBuf[i]      ← serial, cheap

#pragma omp parallel for       ← single fork/join
for pair in index_grm_pairs:
    for i in 0..numNblock:
        N_thread(sampleMissBuf[i])
```

This reduces the OpenMP overhead from $2 \times n_{\text{block}}$ barriers to a single barrier. It also keeps each thread's `N[]` output range hot in L1/L2 across the inner block loop.

### 1.5 Pre-allocated scratch buffers

```cpp
// Constructor:
validIndexBuf.reserve(nMarkerBlock);
sampleMissBuf.resize(numNSampleBlock * markerPerN * numNblockMax);
```

The old code `new uintptr_t[...]` / `delete[]` inside the hot marker loop; the new code reuses pre-allocated buffers, eliminating per-block heap allocation.

### 1.6 `std::popcount` replaces custom wrapper

```cpp
// Before
uint32_t popcounts(uint64_t dw){ return popcount(dw); }  // with __attribute__((target(...)))

// After
sub_miss[k] += std::popcount(block_buf[k]);  // C++20, inlines to POPCNT
```

---

## 2. REML / Heritability Estimation (`main/est_hsq.cpp`)

REML maximises the restricted log-likelihood:

$$\ell(\boldsymbol{\sigma}^2) = -\frac{1}{2}\left(\log|V| + \log|X^T V^{-1} X| + y^T P y\right)$$

where $V = \sum_{i=1}^c \sigma_i^2 A_i + \sigma_e^2 I$ and $P = V^{-1} - V^{-1}X(X^T V^{-1}X)^{-1}X^T V^{-1}$.

Each Newton step requires the score $\partial\ell/\partial\sigma_i^2 = -\frac{1}{2}\bigl(\text{tr}(PA_i) - y^T PA_i Py\bigr)$ and the average-information matrix $H_{ij} = \frac{1}{2}\,y^T PA_i PA_j Py$.

### 2.1 GRM extraction: double loop → Eigen fancy indexing

**Before** — manual O(n²) loop with lower-triangle index arithmetic:
```cpp
#pragma omp parallel for
for (int i = 0; i < _n; i++)
    for (int j = 0; j <= i; j++)
        (_A[0])(j, i) = (_A[0])(i, j) = _grm(kp[i], kp[j]);
```

**After** — symmetrisation then Eigen gather:
```cpp
eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
(_A[0]) = grm_sym(kp_idx, kp_idx);
```

The index into `kp_idx` also removes the `kp[i] >= kp[j]` conditional needed to stay in the lower triangle of the stored GRM.

### 2.2 Diagonal normalisation: loop → vectorised broadcast

Each GRM component $A$ is scaled so $A_{ii} = 1$ for all $i$. This is a column-wise division:

$$A \leftarrow A \cdot D^{-1}, \qquad D = \mathrm{diag}(A)$$

**Before** — O(n²) parallel loop:
```cpp
#pragma omp parallel for
for (int i = 0; i < _n; i++)
    for (int j = 0; j <= i; j++) {
        (_A[0])(i, j) /= (_A[0])(i, i);
        (_A[0])(j, i) = (_A[0])(i, j);
    }
```

**After** — single vectorised column broadcast:
```cpp
eigenVector diag_inv = (_A[0]).diagonal().cwiseInverse();
(_A[0]) = ((_A[0]).array().colwise() / diag_inv.array()).matrix();
(_A[0]).triangularView<Eigen::StrictlyUpper>() =
    (_A[0]).triangularView<Eigen::StrictlyLower>().transpose();
```

### 2.3 V matrix assembly: full matrix → lower triangle only

$V = \sum_i \sigma_i^2 A_i$ is symmetric positive definite. LAPACK `dpotrf` (Cholesky) reads only the lower triangle, so writing the upper half is pure waste: $\frac{n(n-1)}{2}$ needless stores.

**Before** — fills the full n×n matrix:
```cpp
Vi = eigenMatrix::Zero(_n, _n);
for (int i = 0; i < _r_indx.size(); i++)
    Vi += _A[_r_indx[i]] * prev_varcmp[i];
```

**After** — accumulates only the lower triangle (which is all `dpotrf` reads):
```cpp
Vi.triangularView<Eigen::Lower>().setZero();
for (int ci = 0; ci < num_comp; ci++) {
    const double s = prev_varcmp[ci];
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < _n; j++)
        Vi.col(j).tail(_n - j) += s * Ai.col(j).tail(_n - j);
}
```

This halves memory bandwidth for the assembly pass and keeps the columns contiguous in cache.

### 2.4 P matrix: general DGEMM → DSYRK via rank-update

The projection matrix is $P = V^{-1} - V^{-1}X(X^T V^{-1}X)^{-1}X^T V^{-1}$. Let $L = \text{chol}(X^T V^{-1}X)$ (a tiny $c \times c$ Cholesky factor), then with $Z = V^{-1}X L^{-T}$ (an $n \times c$ matrix):

$$P = V^{-1} - Z Z^T$$

This is a symmetric rank-$c$ update — exactly what BLAS DSYRK is designed for.

**Before** — two full n×n DGEMM calls:
```cpp
P = Vi - Vi_X * Xt_Vi_X_i * Vi_X.transpose();
```

**After** — Cholesky of the tiny correction factor, then `DSYRK`:
```cpp
Eigen::LLT<eigenMatrix> llt(Xt_Vi_X_i);   // c×c where c = num_covariates (tiny)
eigenMatrix Z;
Z.noalias() = Vi_X * llt.matrixL();        // n×c
P->triangularView<Eigen::Lower>() = Vi.triangularView<Eigen::Lower>();
P->selfadjointView<Eigen::Lower>().rankUpdate(Z, -1.0);  // DSYRK
P->triangularView<Eigen::Upper>() = P->transpose();
```

`rankUpdate` computes $P_{\text{lower}} = V^{-1}_{\text{lower}} - Z Z^T$, writing only the lower triangle — roughly half the flops and memory bandwidth of the general product.

### 2.5 Factorize-only path (skip-P mode for Hutch++)

When the Hutch++ trace estimator is active, only $\log|V|$ (for the likelihood value) and the ability to solve $Vx = b$ are needed — not the full dense $V^{-1}$. The log-determinant from the Cholesky factor $V = LL^T$ is:

$$\log|V| = 2 \sum_{i=1}^n \log L_{ii}$$

When stochastic trace approximation is active, `calcu_Vi` accepts a `factorize_only = true` flag. The V matrix is Cholesky-factored but `dpotri` is not called, and `Vi` is released immediately:

```cpp
if (factorize_only && !_reml_diagV_adj) {
    Eigen::LLT<eigenMatrix> llt(Vi.selfadjointView<Eigen::Lower>());
    logdet = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
    _Vi_llt = std::move(llt);
    _Vi_use_llt = true;
    Vi.resize(0, 0);   // free n×n memory
    return true;
}
```

This avoids $O(n^3)$ inversion flops and $O(n^2)$ memory for $V^{-1}$ in iterations that only need $\log|V|$.

### 2.6 Implicit P matrix-vector product

Rather than materialising the $n \times n$ P matrix before applying it, an `applyP_vec` function uses cached factors:

$$P v = V^{-1} v - V^{-1} X \underbrace{(X^T V^{-1} X)^{-1}}_{\text{c×c}} X^T V^{-1} v$$

```cpp
eigenVector gcta::applyP_vec(const eigenVector &v) const {
    eigenVector w = _Vi_use_llt
        ? _Vi_llt.solve(v)                                    // triangular solve
        : eigenVector(_Vi.selfadjointView<Eigen::Lower>() * v);  // DSYMV
    eigenVector a = _Vi_X.transpose() * v;    // c-dim GEMV
    eigenVector b = _Xt_Vi_X_i.selfadjointView<Eigen::Lower>() * a;  // c×c SYMV
    w.noalias() -= _Vi_X * b;
    return w;
}
```

Used in AI-REML and Hutch++ to avoid loading the full P matrix from RAM.

### 2.7 Hutch++ stochastic trace estimator

**Motivation**: computing $\text{tr}(PA_i)$ exactly requires either materialising $P$ ($O(n^2)$ storage) or an $O(n^3)$ matrix product. For large $n$ this is the dominant memory bottleneck.

**Standard Girard–Hutchinson estimator**: $\text{tr}(M) \approx \frac{1}{s}\sum_{i=1}^s r_i^T M r_i$ — variance $O(\|M\|_F^2 / s)$.

**Hutch++ estimator** (Meyer et al., 2021): given a rank-$r$ sketch $Q$ of $M$ formed from Rademacher vectors $\{s_i\}$, and $r$ additional probe vectors $\{g_j\}$:

$$\text{tr}(M) \approx \text{tr}(Q^T M Q) + \frac{1}{r}\sum_{j=1}^r g_j^T(I - QQ^T)\,M\,(I - QQ^T)\,g_j$$

Variance decays as $O\!\left(\|M\|_F^2/r^2\right)$ — a quadratically better rate than Girard–Hutchinson.

The probe vectors are generated **once** before the REML loop and reused across all Newton iterations, so the trace estimate is a deterministic function of the variance components — preserving Newton convergence of AI-REML.

**Effect**: eliminates the $8n^2$-byte $P$ matrix from peak memory. Each trace evaluation costs $O(k \cdot n^2)$ (matrix-vector products via DSYMV) instead of $O(n^2)$ storage.

### 2.8 AI matrix `calcu_Hi`: nested loop → Frobenius inner product

The AI matrix entry is $H_{ij} = \frac{1}{2}\,\text{tr}(PA_i\,PA_j)$. Using the Frobenius inner product identity:

$$\text{tr}(BC) = \langle B,\, C^T \rangle_F = \sum_{k,l} B_{kl}\, C_{lk}$$

so $\text{tr}(PA_i\,PA_j) = \langle PA_i,\,(PA_j)^T \rangle_F$, which is the sum of the element-wise (Hadamard) product of $PA_i$ and $(PA_j)^T$.

**Before** — manual nested loop with reduction array:
```cpp
double *d_bufs = new double[_n];
#pragma omp parallel for
for (int k = 0; k < _n; k++)
    for (int l = 0; l < _n; l++)
        d_bufs[k] += PA_i(k, l) * PA_j(l, k);
// ... reduce d_bufs
```

**After** — one Eigen expression:
```cpp
Hi(i, j) = 0.5 * (PA[i].cwiseProduct(PA[j].transpose())).sum();
```

No heap allocation; `cwiseProduct` + `sum` dispatches to a BLAS-level dot-product kernel.

### 2.9 EM-REML: `A*Py` via `selfadjointView`

**Before**:
```cpp
R(i) = (Py.transpose() * _A[_r_indx[i]] * Py)(0, 0);
```
This reads the full symmetric matrix twice.

**After**:
```cpp
tmp.noalias() = _A[_r_indx[i]].selfadjointView<Eigen::Lower>() * Py;  // DSYMV
R(i) = Py.dot(tmp);
```
`selfadjointView` routes to BLAS `DSYMV`, reading only the lower triangle — halving memory bandwidth.

### 2.10 Phenotype covariate adjustment: scalar loop → cast + map

**Before**:
```cpp
eigenVector y_buf = _y.array() - (_X * _b).array();
for (int i = 0; i < n; i++) y[i] = y_buf[i];
```

**After**:
```cpp
Eigen::Map<Eigen::VectorXf>(y.data(), n) = (_y - _X * _b).cast<float>();
```

Single vectorised cast-and-assign, no temporary `double` vector.

---

## 3. MLM Association Tests (`main/mlm_assoc.cpp`)

The mixed-model association test fits, for each candidate SNP $j$:

$$y = x_j \beta_j + C\gamma + g + \varepsilon, \quad g \sim \mathcal{N}(0,\,\sigma_g^2 G),\quad \varepsilon \sim \mathcal{N}(0,\,\sigma_e^2 I)$$

giving the GLS estimator and Wald statistic:

$$\hat{\beta}_j = \frac{x_j^T V^{-1} y}{x_j^T V^{-1} x_j}, \qquad T_j^2 = \frac{\hat{\beta}_j^2}{\widehat{\operatorname{Var}}(\hat{\beta}_j)} \sim \chi^2_1$$

### 3.1 `mlma_calcu_stat`: per-SNP CBLAS → block Eigen GEMM

**Before** — per-SNP matrix-vector products via raw CBLAS:
```cpp
float *Vi  = new float[n*n];
float *X   = new float[n];
float *Vi_X = new float[n];
for (i = 0; i < m; i++) {
    // extract X_snp column
    cblas_sgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, Vi, n, X, 1, 0.0, Vi_X, 1);
    Xt_Vi_X = cblas_sdot(n, X, 1, Vi_X, 1);
    beta[i] = (1.0/Xt_Vi_X) * cblas_sdot(n, y, 1, Vi_X, 1);
}
delete[] Vi; delete[] X; delete[] Vi_X;
```

Per-SNP cost: two SGEMV calls ($2n^2$ flops) + two SDOT calls.

**After** — block GEMM + precomputed `Vi*y`:
```cpp
// Symmetrise Vi once; compute Vi*y once
Eigen::MatrixXf Vi = _Vi.cast<float>().selfadjointView<Eigen::Upper>();
Eigen::VectorXf Vi_y = Vi * y_vec;   // precomputed

// Per block (bs SNPs at once):
Vi_X_block.leftCols(bs).noalias() = Vi * X_block;          // SGEMM: n×bs
xvx_diag.head(bs) = (X_block.cwiseProduct(Vi_X_block.leftCols(bs))).colwise().sum(); // diagonal only
Xt_Vi_y_block.head(bs).noalias() = X_block.transpose() * Vi_y;  // SGEMV
```

Key improvements:
- `Vi * y_vec` hoisted out of the m-SNP loop (was recomputed per SNP via `sdot`)
- Per-block SGEMM instead of m per-SNP SGEMV calls, giving much better BLAS utilisation
- `diag(X^T Vi X)` computed via elementwise product + column sum — only the diagonal is needed, avoiding an O(bs²) matrix multiply

### 3.2 `mlma_calcu_stat_covar`: Schur complement precomputation

The covariate case previously solved a $(p+1) \times (p+1)$ augmented system per SNP by rebuilding `Xt_Vi_X` each time, costing $O(p^3)$ per SNP.

**New approach** — partition the covariate block using the Schur complement. Let $A = C^T V^{-1} C$ (a tiny $p \times p$ matrix). For each SNP $j$, define $d_j = C^T V^{-1} x_j$. The Schur complement of $A$ in the full $(p+1)\times(p+1)$ system is:

$$S_j = x_j^T V^{-1} x_j - d_j^T A^{-1} d_j$$

which equals $1/\operatorname{Var}(\hat{\beta}_j)$. The effect estimate is:

$$\hat{\beta}_j = \frac{f_j - d_j^T A^{-1} t}{S_j}, \qquad f_j = x_j^T V^{-1} y, \quad t = A^{-1} C^T V^{-1} y$$

**Implementation** — $A$, $A^{-1}$ (via LLT), and $t$ are **precomputed once** before the SNP loop; only scalar operations on precomputed vectors remain per block:

1. **Precompute once**: $V_C = V^{-1}C$ (n×p GEMM), $A = C^T V_C$ (p×p), $\text{LLT}(A)$, $t = A^{-1}(C^T V^{-1} y)$
2. **Per block of $b_s$ SNPs**:
   - $V_X = V^{-1} X_{\text{block}}$ (n×$b_s$ GEMM)
   - $D = C^T V_X$ (p×$b_s$)
   - $E = A^{-1} D$ (p×$b_s$ triangular solve)
   - Diagonal Schur: $S_j = \text{diag}(X^T V_X) - \text{diag}(D^T E)$ — avoids any $O(b_s^2)$ product

The per-SNP cost drops from an $O(p^3)$ matrix inversion + two O(n²) GEMMs to a handful of dot products on pre-computed quantities.

```cpp
// Precompute (once, before SNP loop)
Eigen::MatrixXf Vi_C(n, p);  Vi_C.noalias() = Vi * C;
A_mat.noalias() = Vi_C.transpose() * C;   // p×p
Eigen::LLT<Eigen::MatrixXf> A_llt(A_mat);
t_vec = A_llt.solve(C.transpose() * Vi_y);
Vi_C.resize(0, 0);  // free n×p memory before the loop

// Per block
Vi_X_block.leftCols(bs).noalias() = Vi * X_block;
D_block.leftCols(bs).noalias() = C.transpose() * Vi_X_block.leftCols(bs);
E_block.leftCols(bs).noalias() = A_llt.solve(D_block.leftCols(bs));
// Diagonal Schur complement: avoid O(bs²) full matrix product
xvx_diag.head(bs)     = (X_block.cwiseProduct(Vi_X_block.leftCols(bs))).colwise().sum();
d_dot_e_diag.head(bs) = (D_block.leftCols(bs).cwiseProduct(E_block.leftCols(bs))).colwise().sum();
// S = diag(X^T Vi X) - diag(D^T E)
```

### 3.3 Vi symmetrisation

In both stat functions, `_Vi` (stored as lower-triangular) is symmetrised into a plain dense matrix exactly once before the SNP loop:

```cpp
Eigen::MatrixXf Vi = _Vi.cast<float>().selfadjointView<Eigen::Upper>();
_Vi.resize(0, 0);  // release double-precision copy immediately
```

This ensures all downstream GEMMs dispatch through Eigen's full BLAS path. `selfadjointView` expressions as the *target* of a block product can fall back to a scalar kernel in some Eigen/BLAS configurations; a plain dense matrix never does.

---

## 4. Further Algorithmic Improvements

### 4.1 LD pruning: sequential dot products → block GEMM (`LinAlg.cpp`)

Pairwise squared correlations within a block of $b$ SNPs are:

$$R^2_{ij} = \left(\frac{z_i^T z_j}{n}\right)^2 = \left[\frac{1}{n} Z^T Z\right]_{ij}^2$$

where $z_i$ are standardised genotype vectors. The full block $R^2$ matrix is the Hadamard (element-wise) square of $\frac{1}{n}Z^T Z$.

**Before** — nested loop computing each $r^2_{ij}$ sequentially as a dot product.

**After** — single DGEMM then element-wise square:
```cpp
eigenMatrix R = (X_sub.transpose() * X_sub) * (1.0 / n);  // DGEMM: b×b
R = R.array().square();                                     // r²_ij = R_ij²
```
One level-3 BLAS call replaces $O(b^2)$ sequential level-1 calls.

### 4.2 Approximate PCA: Spectra Lanczos / randomised SVD (`grm.cpp`)

Full `SelfAdjointEigenSolver` costs $O(n^3)$ flops. Two approximate paths were added:

**Lanczos (`--pca-approx Lanczos`)** — implicitly-restarted Lanczos via Spectra; $O(k \cdot n^2)$ total:

$$\frac{\text{cost}_{\text{Lanczos}}}{\text{cost}_{\text{full}}} \approx \frac{k}{n}$$

For $n = 10^5$, $k = 20$: ~5000× fewer flops.

```cpp
Spectra::DenseSymMatProd<double> op(G);
int ncv = std::min(n, 2 * k + 10);
Spectra::SymEigsSolver<...> eigs(op, k, ncv);
eigs.compute(Spectra::SortRule::LargestAlge);
```

**Randomised SVD (`--pca-approx SVD`)** — random sketch with power iteration. With $q$ power iterations, error decays as $O(\sigma_{k+1}^{2q+1})$:
```cpp
for (int i = 0; i < pca_power_iter; ++i)
    sketch = G * (G * sketch);   // power iteration
```

### 4.3 Median: $O(n \log n)$ → $O(n)$ selection (`CommFunc.cpp`)

**Before** — `std::stable_sort`. **After** — `std::nth_element` (O(n) average). Used in LD analysis reporting median $r^2$ across thousands of genomic windows.

### 4.4 Log-determinant from factor diagonal (`Matrix.hpp`, `LinAlg.cpp`)

Previously called `lu.determinant()` (can underflow/overflow for large matrices) and then `lu.logAbsDeterminant()` — two decompositions. Now extracted in a single pass:

$$\log|V| = 2\sum_{i=1}^n \log L_{ii} \quad (LLT), \qquad \log|\det A| = \sum_i \log|U_{ii}| \quad (LU)$$

```cpp
// LLT
logdet = llt.matrixL().nestedExpression().diagonal().array().log().sum() * 2.0;
// LU
logdet = lu.matrixLU().diagonal().array().abs().log().sum();
```

---

## Performance Summary Table

| File | Change | Technique | Impact |
|------|--------|-----------|--------|
| GRM.cpp | DSYRK via `rankUpdate` | Symmetry exploit | Write $n(n+1)/2$ vs $n^2$ elements |
| GRM.cpp | `stdGenoLD` 64-byte alignment | Memory layout | Enables aligned AVX-512 BLAS kernel |
| GRM.cpp | `nMarkerBlock` 128 → 1024 | Tuning | 8× fewer rank-$k$ dispatch calls |
| GRM.cpp | OpenMP hoist: $2n_b$ → 1 barrier | Loop hoisting | Eliminates repeated fork/join overhead |
| GRM.cpp | Pre-alloc `sampleMissBuf` | Memory | No per-block heap allocation |
| GRM.cpp | `std::popcount` | Instruction | Inlines to single `POPCNT` |
| est_hsq | GRM gather via fancy index | Vectorisation | Removes scalar `j ≤ i` conditional |
| est_hsq | Diagonal normalise: `colwise` | Vectorisation | $A \leftarrow A D^{-1}$ in one BLAS pass |
| est_hsq | V assembly: lower triangle only | Symmetry exploit | ~$2\times$ fewer writes |
| est_hsq | P via `rankUpdate` (DSYRK) | Symmetry exploit | ~$2\times$ flops + bandwidth |
| est_hsq | Factorize-only: skip $V^{-1}$ | Algorithm | Saves $O(n^3) + O(n^2)$ in Hutch++ mode |
| est_hsq | `applyP_vec`: implicit $Pv$ | Algorithm | $O(cn)$ per call vs $O(n^2)$ load |
| est_hsq | Hutch++ $\text{tr}(PA_i)$ | Algorithm | Eliminates $8n^2$-byte P matrix from memory |
| est_hsq | AI matrix: Frobenius inner product | Vectorisation | Replaces O(n²) loop; no heap alloc |
| est_hsq | EM-REML: DSYMV via `selfadjointView` | Symmetry exploit | ~$2\times$ bandwidth for $A_i v$ |
| mlm_assoc | $V^{-1}y$ hoisted outside SNP loop | Loop hoisting | Saves $m \times 2n^2$ flops |
| mlm_assoc | Per-SNP SGEMV → block SGEMM | Cache reuse | Level-3 vs level-2 BLAS |
| mlm_assoc | Diagonal-only $x^T V^{-1} x$ | Avoid work | $O(n b_s)$ vs $O(n b_s + b_s^2)$ |
| mlm_assoc | Schur complement precomputed once | Loop hoisting | $O(p^3)$ per SNP → $O(p b_s)$ per block |
| mlm_assoc | $V^{-1}$ freed after symmetrisation | Memory | Releases $8n^2$ bytes; guarantees BLAS fast path |
| LinAlg.cpp | LD block GEMM $r^2$ | Level-3 BLAS | $O(b^2)$ dot products → one DGEMM |
| grm.cpp | Spectra Lanczos PCA | Algorithm | $O(kn^2)$ vs $O(n^3)$ eigensolver |
| grm.cpp | Randomised SVD + power iteration | Algorithm | $O(kn^2)$; better spectral accuracy |
| CommFunc.cpp | `median`: `nth_element` | Algorithm | $O(n \log n)$ → $O(n)$ |
| Matrix.hpp | Logdet from factor diagonal | Numerical | One decomposition; no overflow |

---

# Part II: Supporting Changes

These changes improve portability, correctness, and I/O. They are not individually performance changes, but several remove the constraints (x86-only build, MKL requirement, OMP spin-polling) that previously prevented the Part I improvements from delivering their full benefit.

---

## 5. MKL Removal and BLAS Portability

**`main/mkl.cpp` deleted** (~800 lines). It contained `make_grm_mkl`, `make_XMat_mkl`, `std_XMat_mkl`, `std_XMat_d_mkl` — all wrapped in `#if GCTA_CPU_x86` guards, making the GRM pipeline x86-Linux-only. The Eigen-based paths in §§1–3 now serve all platforms (x86, Apple Silicon, ARM Linux) with any BLAS backend.

All remaining `#if GCTA_CPU_x86` guards for LAPACK/BLAS ABI selection (Fortran `foo_` vs C `foo`) were also removed — Eigen handles the ABI internally.

**BLAS backend selection** (`CMakeLists.txt`):
```cmake
set(GCTA_BLAS_BACKEND "Auto" CACHE STRING
    "BLAS/LAPACK backend: Auto, OpenBLAS, MKL, or Accelerate")
```

**Native optimisation** — `-march=native` (x86) / `-mcpu=native` (Apple Silicon) enabled by default in Release builds, plus IPO/LTO. This makes the aligned `stdGenoLD` stride and Eigen expression templates immediately useful on every machine.

**OMP spin-wait elimination** (`main.cpp`):
```cpp
setenv("KMP_BLOCKTIME", "0", 0);   // Intel OMP: sleep after parallel region
setenv("GOMP_SPINCOUNT", "0", 0);  // GNU OMP: yield after parallel region
```
Without these, OMP threads spin-poll after each `#pragma omp parallel` region, consuming memory bandwidth and starving the immediately following BLAS call.

---

## 6. Statistical Libraries: Boost.Math

**~100 lines of hand-rolled numerical series** (`betai`, `betacf`, `gammln`, `gammp`, `gser`, `gcf`) deleted from `StatFunc.cpp` and replaced with Boost.Math:
```cpp
boost::math::cdf(boost::math::complement(boost::math::chi_squared(df), x))
```
Boost.Math uses Lanczos-approximated special functions with guaranteed relative error bounds.

A **log-space path** was added for extreme $\chi^2$ values (`--log-pval`): reports $\ln(p)$ using the asymptotic expansion $\ln Q(a,x) \approx -x + (a-1)\ln x - \ln\Gamma(a)$ for $x \gg a$, avoiding underflow to $p = 0$ for GWAS hits below $10^{-300}$.

`StatLib.cpp`'s **`rankContrast`** replaced 40 lines of `dgeqrf`/`dormqr` LAPACK calls guarded by `#if GCTA_CPU_x86` with 8 lines of `Eigen::HouseholderQR`.

---

## 7. I/O and Memory Management

| Change | Location | Before | After |
|--------|----------|--------|-------|
| gzip read/write | `grm.cpp` | Custom `gzifstream`/`gzofstream` (~200-line C wrapper) | `boost::iostreams::filtering_stream` |
| FastFAM covar conditioning | `FastFAM.cpp` | 2 CBLAS `dgemv` + 10-line `#if GCTA_CPU_x86` | `y -= covar * (H * y)` (1 line) |
| GRM member `sub_miss` | `GRM.h` / `GRM.cpp` | `uint32_t*` raw pointer, manual `delete[]` | `std::vector<uint32_t>` |
| FastFAM scratch arrays | `FastFAM.cpp` | `new double[n]` × 4, manual `delete[]` | `std::vector<double>` × 4 |
| `quantile` | `CommFunc.cpp` | Mapped `vector` into `Eigen::VectorXd` | Direct `quantile(data, n, prob)` — no temporary |
| `rand_seed()` | `CommFunc.cpp` | Reversed time string, `abs(atoi(...))` | `std::random_device{}() & 0x7FFFFFFF` |

---

## 8. New CLI Flags (Activating Part I Features)

| Flag | Activates | Effect |
|------|-----------|--------|
| `--reml-trace-approx` | Hutch++ (§2.7) | Stochastic $\text{tr}(PA_i)$; eliminates $O(n^2)$ P matrix |
| `--reml-trace-nprobes N` | Hutch++ (§2.7) | Number of probe vectors (default 90) |
| `--save-reml FILE` / `--load-reml FILE` | Checkpoint | Skip $O(n^3)$ REML re-fit in MLMA-LOCO |
| `--log-pval` | Log-space $p$ (§6) | Reports $\ln(p)$; avoids underflow for $p < 10^{-300}$ |
| `--pca-approx Lanczos\|SVD` | Spectra / rSVD (§4.2) | $O(kn^2)$ vs $O(n^3)$ PCA for large GRMs |

---

## References

- Meyer, R. A., Musco, C., Musco, C., & Woodruff, D. P. (2021). *Hutch++: Optimal Stochastic Trace Estimation*. SIAM Symposium on Simplicity in Algorithms (SOSA).
- Halko, N., Martinsson, P.-G., & Tropp, J. A. (2011). *Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions*. SIAM Review.
- Eigen documentation: `selfadjointView`, `rankUpdate`, `triangularView`, `LLT`, `PartialPivLU`.
- Boost.Math documentation: *Statistical distributions and special functions*.
- Spectra documentation: `SymEigsSolver` (implicitly-restarted Lanczos).
