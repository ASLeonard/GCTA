/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * RemlCtx — flat workspace struct that carries all state for a single-GRM
 * univariate REML computation.  This mirrors the REML-relevant private
 * members of the legacy `gcta` class but is usable from any v2 translation
 * unit without instantiating a `gcta` object.
 *
 * Uses Eigen types directly (double precision throughout) so it can be
 * included without pulling in gcta.h.
 *
 * Scope: single-GRM, univariate, non-bivariate, non-within-family.
 * Bivariate / within-family REML remain in gcta and are out of scope here.
 */
#pragma once

#include "mlma_woodbury.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

// REML is always double-precision irrespective of SINGLE_PRECISION build flag.
using RemlMat = Eigen::MatrixXd;
using RemlVec = Eigen::VectorXd;

// ──────────────────────────────────────────────────────────────────────────────
// RemlCtx: input + config + workspace + output
// Caller fills the "Input" and "Config" sections, then calls reml::compute().
// After compute() returns, "Output" holds variance components and fixed effects;
// "Workspace" holds the V^{-1} state needed for the subsequent association test.
// ──────────────────────────────────────────────────────────────────────────────
struct RemlCtx {
    // ── Input (set by caller before reml::compute()) ──────────────────────────
    int  n     = 0;   // sample size
    int  X_c   = 0;   // columns in design matrix (1 + n_covariates)
    RemlMat  X;       // n × X_c (intercept in column 0)
    RemlVec  y;       // phenotype vector (n)
    double y_Ssq = 0.0; // sample variance of y (computed by caller before call)

    // GRM components: A[r_indx[i]] is the i-th variance component matrix.
    // For standard single-GRM REML: A = {GRM, empty} and r_indx = {0, 1}.
    // An empty matrix (size() == 0) represents the identity component (residual).
    std::vector<RemlMat> A;
    std::vector<int>     r_indx;  // typically {0, 1}
    RemlMat grm_N;   // from _grm_N — used only for Woodbury auto-k M estimation

    // ── Config (algorithm control) ────────────────────────────────────────────
    int    reml_mtd                  = 0;     // 0=AI-REML, 1=Fisher, 2=EM-REML
    int    reml_max_iter             = 100;
    int    reml_inv_mtd              = 0;     // 0=LLT, 1=LU (for V inversion fallback)
    double reml_diag_mul             = 0.01;  // diagonal jitter when V is near-singular
    int    reml_diagV_adj            = 0;     // 0=none, 1=diag jitter, 2=bending
    int    woodbury_rank             = 0;     // >0 fixed k; -1 auto-k; 0 off
    double woodbury_buffer_factor    = 1.5;   // MP signal buffer for auto-k
    int    woodbury_k_max            = 0;     // rank cap for auto-k (0 → min(n−1,1200))
    bool   woodbury_nystrom          = false; // true → single-pass Nystrom basis
    bool   reml_trace_approx         = false; // Hutch++ trace (skips n x n P)
    int    reml_trace_approx_nprobes = 90;
    int    reml_trace_power_iter     = 0;     // power-iter for Hutch++ range sketch
    bool   reml_force_inv            = false;
    bool   reml_force_converge       = false;
    bool   reml_no_converge          = false;
    bool   reml_fixed_var            = false;
    bool   reml_allow_constrain_run  = false;

    // Output naming (for .hsq file and LOGGER lines)
    std::string out;                       // output file prefix
    std::vector<std::string> var_name;     // e.g. {"V(G)", "V(e)"}
    std::vector<std::string> hsq_name;     // e.g. {"V(G)/Vp"}

    // ── Woodbury basis (set by reml::compute_woodbury_basis()) ───────────────
    bool    Vi_use_woodbury = false;
    int     woodbury_rank_  = 0;     // actual rank used (<= woodbury_rank or auto)
    RemlMat Uk;                      // n x k leading eigenvectors of K
    RemlVec dk;                      // k eigenvalues (clamped >= 0)
    double  lambda_tail     = 0.0;   // average bulk eigenvalue
    double  tail_d_var      = 0.0;   // var of tail eigenvalues (2nd-order logdet corr)
    double  sigma2_eff      = 0.0;   // sigma2_g * lambda_tail + sigma2_e (updated per AI iter)
    double  sg2             = 0.0;   // sigma2_g (cached for calcu_tr_PA_woodbury)
    RemlVec ck;                      // Woodbury corrections c_j

    // ── Iteration workspace (managed by reml::compute()) ─────────────────────
    bool    Vi_use_llt = false; // true when Vi_L is valid (factorize_only path)
    RemlMat Vi;                 // full V^{-1} (or empty when Woodbury)
    RemlMat Vi_L;               // Cholesky L from in-place dpotrf
    RemlMat Vi_X;               // V^{-1} X  (n x X_c)
    RemlMat Xt_Vi_X_i;          // (X' V^{-1} X)^{-1}  (X_c x X_c)
    RemlMat Uk_Vi_X;            // U_k^T V^{-1} X  (k x X_c, Woodbury)
    RemlMat UkTX;               // U_k^T X  (k x X_c, constant, cached)
    RemlVec UkTy;               // U_k^T y  (k, constant, cached)
    RemlMat hutchpp_S;          // Hutch++ Rademacher probes  (n x k)
    RemlMat hutchpp_G;          // Hutch++ Rademacher probes  (n x k)
    RemlMat P;                  // projection matrix (freed after convergence)

    // ── Status flags ─────────────────────────────────────────────────────────
    bool reml_AI_not_invertible = false;
    bool reml_have_bend_A       = false;

    // ── Output (populated by reml::compute()) ────────────────────────────────
    std::vector<double> varcmp; // variance components in r_indx order
    RemlVec b;                  // fixed-effect estimates (X_c elements)

    // Minimal REML summary for MLMA-style .hsq output.
    // Mirrors the non-bivariate mlmassoc summary fields in main/est_hsq.cpp.
    std::vector<double> varcmp_se; // SE for each variance component (sqrt(diag(Hi)))
    double Vp = 0.0;               // total phenotypic variance
    double Vp_se = 0.0;            // SE(Vp)
    std::vector<double> hsq;       // V(G_i)/Vp for non-residual components
    std::vector<double> hsq_se;    // SE for hsq entries
    double logL = 0.0;             // converged REML log-likelihood
    bool has_logL = false;
};
