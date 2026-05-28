#pragma once
/*
 * Shared Woodbury-MLMA types and helpers used by both the classic gcta MLMA
 * path (main/mlm_assoc.cpp) and the streaming MLMA path (src/MLMA_stream.cpp).
 *
 * Layout convention: Uk_f is stored as k×n (rows=k, cols=n) — transposed
 * relative to the n×k eigenvector matrix — so that Uk_f * x  produces a k×1
 * product without an explicit .transpose() call, which is cache-friendlier for
 * the n×BLOCK batch products in the hot MLMA loop.
 *
 * All helper functions below assume this k×n layout.
 */

#include <Eigen/Dense>
#include <cmath>
#include "main/StatFunc.h"

// ---------------------------------------------------------------------------
// WoodburyMLMACache
// Float-precision factors sufficient to apply V^{-1} via the Woodbury identity.
// Uk_f is k×n (transposed for GEMM efficiency).
// ---------------------------------------------------------------------------
struct WoodburyMLMACache {
    Eigen::MatrixXf Uk_f;          // k×n  (transposed eigenvectors)
    Eigen::VectorXf ck_f;          // k    (Woodbury coefficients)
    Eigen::VectorXf sqrt_ck_f;     // k    cwiseSqrt(ck_f); safe because ck >= 0
    float           sigma2_eff_f = 0.0f;
};

// ---------------------------------------------------------------------------
// woodbury_apply_Vi_f
// Apply V^{-1} to a float vector:
//   V^{-1} v = (v − Uk^T (ck ⊙ (Uk v))) / σ²_eff
// Uk is k×n so Uk*v produces k×1 directly (no explicit transpose needed).
// ---------------------------------------------------------------------------
inline Eigen::VectorXf woodbury_apply_Vi_f(const WoodburyMLMACache& wb,
                                            const Eigen::VectorXf& v)
{
    Eigen::VectorXf Ukv = wb.Uk_f * v;                    // k×n * n×1 = k×1
    Ukv.array() *= wb.ck_f.array();
    return (v - wb.Uk_f.transpose() * Ukv) / wb.sigma2_eff_f;  // n×1
}

// ---------------------------------------------------------------------------
// woodbury_apply_Vi_mat_f
// Apply V^{-1} to each column of a float matrix.
// ---------------------------------------------------------------------------
inline Eigen::MatrixXf woodbury_apply_Vi_mat_f(const WoodburyMLMACache& wb,
                                                const Eigen::MatrixXf& M)
{
    Eigen::MatrixXf UkM = wb.Uk_f * M;                    // k×n * n×p = k×p
    UkM.array().colwise() *= wb.ck_f.array();
    return (M - wb.Uk_f.transpose() * UkM) / wb.sigma2_eff_f;  // n×p
}

// ---------------------------------------------------------------------------
// woodbury_xvx_diag_block
// Compute diag(X^T V^{-1} X) for the first bs columns of X_block using:
//   x^T V^{-1} x = (||x||² − ||sqrt_ck ⊙ (Uk x)||²) / σ²_eff
// Result is clamped to >= 0 to guard against floating-point cancellation.
// ---------------------------------------------------------------------------
inline void woodbury_xvx_diag_block(const WoodburyMLMACache& wb,
                                     const Eigen::MatrixXf& X_block, Eigen::Index bs,
                                     Eigen::VectorXf& xvx_diag)
{
    Eigen::MatrixXf UkX = wb.Uk_f * X_block.leftCols(bs);  // k×n * n×bs = k×bs
    UkX.array().colwise() *= wb.sqrt_ck_f.array();
    xvx_diag.head(bs) = (X_block.leftCols(bs).colwise().squaredNorm()
                         - UkX.colwise().squaredNorm()) / wb.sigma2_eff_f;
    xvx_diag.head(bs) = xvx_diag.head(bs).cwiseMax(0.0f);
}

// ---------------------------------------------------------------------------
// mlma_snp_stat
// Compute beta, se, and p-value for a single SNP.
//
//   numerator   = x^T Vi y  (or x^T Vi y − covariate correction in the
//                             covariate-adjusted path)
//   denominator = x^T Vi x  (or its Schur complement in the covariate path)
//
// Returns false if denominator <= 1e-30 (degenerate/monomorphic SNP);
// beta, se, and pval are unchanged in that case.
// ---------------------------------------------------------------------------
[[nodiscard]] inline bool mlma_snp_stat(float numerator, float denominator,
                                         bool log_pval,
                                         float& beta, float& se, double& pval)
{
    if (denominator <= 1.0e-30f) return false;
    const float inv_d    = 1.0f / denominator;
    const float sqrt_inv = std::sqrt(inv_d);
    beta = numerator * inv_d;
    se   = sqrt_inv;
    // chisq = beta/se = (numerator/denom) / (1/sqrt(denom)) = numerator/sqrt(denom)
    const float chisq = numerator * sqrt_inv;
    pval = StatFunc::pchisq(static_cast<double>(chisq * chisq), 1, log_pval);
    return true;
}
