/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * RemlEngine — free-function REML engine for the v2 MLMALoco path.
 * Single-GRM, univariate, non-bivariate, non-within-family REML only.
 *
 * Extracted from main/est_hsq.cpp with the following substitutions:
 *   this->_foo       → ctx.foo
 *   eigenMatrix      → RemlMat  (= Eigen::MatrixXd)
 *   eigenVector      → RemlVec  (= Eigen::VectorXd)
 *   _bivar_reml      → false (branches removed)
 *   _within_family   → false (branches removed)
 */

#include "RemlEngine.hpp"
#include "RemlCtx.hpp"
#include "RemlState.hpp"
#include "Matrix.hpp"
#include "Logger.h"
#include "cpu.h"
#include "mlma_woodbury.hpp"

#include <Eigen/Dense>
#include <random>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#ifdef __APPLE__
  #include <omp.h>
#else
  #include <omp.h>
#endif

// ──────────────────────────────────────────────────────────────────────────────
// Local type aliases (double-precision throughout)
// ──────────────────────────────────────────────────────────────────────────────
using RemlMat = Eigen::MatrixXd;
using RemlVec = Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers (file-local, not exposed via RemlEngine.hpp)
// ─────────────────────────────────────────────────────────────────────────────

namespace {

// Woodbury helpers — mirrors gcta::woodbury_Kv / woodbury_Viv / woodbury_ViZ
// but takes ctx by const ref instead of reading this->_ members.

RemlVec woodbury_Kv(const RemlCtx& ctx, const RemlVec& v) {
    RemlVec Ukv = ctx.Uk.transpose() * v;
    Ukv.array() *= (ctx.dk.array() - ctx.lambda_tail);
    return ctx.lambda_tail * v + ctx.Uk * Ukv;
}

RemlVec woodbury_Viv(const RemlCtx& ctx, const RemlVec& v) {
    RemlVec Ukv = ctx.Uk.transpose() * v;
    Ukv.array() *= ctx.ck.array();
    return (v - ctx.Uk * Ukv) / ctx.sigma2_eff;
}

RemlMat woodbury_ViZ(const RemlCtx& ctx, const RemlMat& Z) {
    RemlMat UkZ = ctx.Uk.transpose() * Z;
    UkZ.array().colwise() *= ctx.ck.array();
    return (Z - ctx.Uk * UkZ) / ctx.sigma2_eff;
}

// Implicit P*v without materialising n×n P
RemlVec applyP_vec(const RemlCtx& ctx, const RemlVec& v) {
    RemlVec w;
    if (ctx.Vi_use_woodbury) {
        w = woodbury_Viv(ctx, v);
    } else if (ctx.Vi_use_llt) {
        w = v;
        ctx.Vi_L.triangularView<Eigen::Lower>().solveInPlace(w);
        ctx.Vi_L.triangularView<Eigen::Lower>().adjoint().solveInPlace(w);
    } else {
        w = RemlVec(ctx.Vi.selfadjointView<Eigen::Lower>() * v);
    }
    RemlVec a = ctx.Vi_X.transpose() * v;
    RemlVec b = ctx.Xt_Vi_X_i.selfadjointView<Eigen::Lower>() * a;
    w.noalias() -= ctx.Vi_X * b;
    return w;
}

RemlMat applyP_mat(const RemlCtx& ctx, const RemlMat& Z) {
    RemlMat W;
    if (ctx.Vi_use_woodbury) {
        W = woodbury_ViZ(ctx, Z);
    } else if (ctx.Vi_use_llt) {
        W = Z;
        ctx.Vi_L.triangularView<Eigen::Lower>().solveInPlace(W);
        ctx.Vi_L.triangularView<Eigen::Lower>().adjoint().solveInPlace(W);
    } else {
        W = RemlMat(ctx.Vi.selfadjointView<Eigen::Lower>() * Z);
    }
    const RemlMat A = ctx.Vi_X.transpose() * Z;
    W.noalias() -= ctx.Vi_X * (ctx.Xt_Vi_X_i.selfadjointView<Eigen::Lower>() * A);
    return W;
}

bool inverse_H(RemlCtx& ctx, RemlMat& H) {
    double d_buf = 0.0;
    INVmethod method = (ctx.reml_inv_mtd == 0) ? INV_LLT : static_cast<INVmethod>(ctx.reml_inv_mtd);
    int rank = 0;
    return SquareMatrixInverse(H, d_buf, rank, method);
}

void bend_V(RemlCtx& ctx, RemlMat& Vi) {
    Eigen::SelfAdjointEigenSolver<RemlMat> eigensolver(Vi);
    RemlVec eval = eigensolver.eigenvalues();
    // bending_eigenval inline:
    double eval_m = eval.mean();
    if (eval.minCoeff() > 1e-6) {
        // no bending needed — just reconstruct
    } else {
        double S = 0.0, P = 0.0;
        for (int j = 0; j < eval.size(); j++) {
            if (eval[j] >= 0) continue;
            S += eval[j];
            P = -eval[j];
        }
        double W = S * S / ctx.reml_diag_mul + 1.0;
        for (int j = 0; j < eval.size(); j++) {
            if (eval[j] >= 0) continue;
            eval[j] = P * (S - eval[j]) * (S - eval[j]) / W;
        }
        eval *= eval_m / eval.mean();
    }
    eval.array() = 1.0 / eval.array();
    Vi = eigensolver.eigenvectors() * eval.asDiagonal() * eigensolver.eigenvectors().transpose();
}

int constrain_varcmp(const RemlCtx& ctx, RemlVec& varcmp) {
    const int m = static_cast<int>(ctx.r_indx.size());
    double delta = 0.0;
    constexpr double constr_scale = 1e-6;
    int num = 0;
    std::vector<int> constrain(m, 0);
    for (int i = 0; i < m; i++) {
        if (varcmp[i] < 0) {
            delta += ctx.y_Ssq * constr_scale - varcmp[i];
            varcmp[i] = ctx.y_Ssq * constr_scale;
            constrain[i] = 1;
            num++;
        }
    }
    if (num < m) {
        delta /= (m - num);
        for (int i = 0; i < m; i++)
            if (constrain[i] < 1 && varcmp[i] > delta) varcmp[i] -= delta;
    }
    return num;
}

void init_varcomp(const RemlCtx& ctx,
                  const std::vector<double>& priors_var,
                  const std::vector<double>& priors,
                  RemlVec& varcmp) {
    const int m = static_cast<int>(ctx.r_indx.size());
    varcmp = RemlVec::Zero(m);

    if (!priors_var.empty()) {
        for (int i = 0; i < m - 1; i++) varcmp[i] = priors_var[i];
        if ((int)priors_var.size() < m)
            varcmp[m - 1] = ctx.y_Ssq - varcmp.sum();
        else
            varcmp[m - 1] = priors_var[m - 1];
    } else if (!priors.empty()) {
        double d_buf = 0.0;
        for (int i = 0; i < m - 1; i++) {
            varcmp[i] = priors[i] * ctx.y_Ssq;
            d_buf += priors[i];
        }
        if (d_buf > 1.0) LOGGER.e(0, "\n  --reml-priors. The sum of all prior values should not exceed 1.0.");
        varcmp[m - 1] = (1.0 - d_buf) * ctx.y_Ssq;
    } else {
        varcmp.setConstant(ctx.y_Ssq / m);
    }

    // Woodbury HE warm-start (single-GRM, no priors)
    if (ctx.Vi_use_woodbury && (int)ctx.r_indx.size() == 2
            && priors.empty() && priors_var.empty()) {
        const int n = ctx.n;
        const int k = ctx.woodbury_rank_;
        const double Vy = ctx.y_Ssq;
        const double trK  = ctx.dk.sum() + static_cast<double>(n - k) * ctx.lambda_tail;
        const double trK2 = ctx.dk.squaredNorm()
                          + static_cast<double>(n - k) * (ctx.tail_d_var + ctx.lambda_tail * ctx.lambda_tail);
        const RemlVec Ky = woodbury_Kv(ctx, ctx.y);
        const double yKy = ctx.y.dot(Ky);
        const double yy  = ctx.y.squaredNorm();
        const double denom = static_cast<double>(n) * trK2 - trK * trK;
        double sg_he = (denom > 1e-10)
            ? (static_cast<double>(n) * yKy - trK * yy) / denom
            : Vy * 0.5;
        sg_he = std::max(sg_he, 0.01 * Vy);
        sg_he = std::min(sg_he, 0.99 * Vy);
        varcmp(0) = sg_he;
        varcmp(1) = std::max(Vy - sg_he, 0.01 * Vy);
    }
}

// Returns false if V is not positive-definite and inversion failed.
bool calcu_Vi(RemlCtx& ctx, RemlVec& prev_varcmp, double& logdet, int& iter, bool factorize_only) {
    ctx.Vi_use_llt = false;

    if (ctx.Vi_use_woodbury) {
        if (!factorize_only)
            LOGGER.e(0, "Woodbury REML is incompatible with explicit V^{-1} materialisation.");
        ctx.Vi.resize(0, 0);
        const double sg2 = prev_varcmp[0];
        const double se2 = prev_varcmp[ctx.r_indx.size() - 1];
        ctx.sigma2_eff = sg2 * ctx.lambda_tail + se2;
        if (ctx.sigma2_eff <= 0.0) return false;
        ctx.sg2 = sg2;
        const int k = ctx.woodbury_rank_;
        ctx.ck.resize(k);
        logdet = static_cast<double>(ctx.n - k) * std::log(ctx.sigma2_eff);
        for (int j = 0; j < k; ++j) {
            const double delta    = std::max(0.0, ctx.dk[j] - ctx.lambda_tail);
            const double sig_delta = sg2 * delta;
            ctx.ck[j] = sig_delta / (ctx.sigma2_eff + sig_delta);
            logdet += std::log(ctx.sigma2_eff + sig_delta);
        }
        if (ctx.tail_d_var > 0.0) {
            const double r = sg2 / ctx.sigma2_eff;
            logdet -= 0.5 * r * r * static_cast<double>(ctx.n - k) * ctx.tail_d_var;
        }
        return true;
    }

    // Dense path: assemble V (lower triangle)
    if (factorize_only && static_cast<int>(ctx.r_indx.size()) > 1)
        ctx.Vi.swap(ctx.Vi_L);
    ctx.Vi.resize(ctx.n, ctx.n);

    if (ctx.r_indx.size() == 1) {
        ctx.Vi.diagonal() = RemlVec::Constant(ctx.n, 1.0 / prev_varcmp[0]);
        logdet = ctx.n * std::log(prev_varcmp[0]);
    } else {
        const int num_comp = static_cast<int>(ctx.r_indx.size());
        ctx.Vi.triangularView<Eigen::Lower>().setZero();
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < ctx.n; j++) {
            for (int ci = 0; ci < num_comp; ci++)
                if (ctx.A[ctx.r_indx[ci]].size() > 0)
                    ctx.Vi.col(j).tail(ctx.n - j) +=
                        prev_varcmp[ci] * ctx.A[ctx.r_indx[ci]].col(j).tail(ctx.n - j);
        }
        for (int ci = 0; ci < num_comp; ci++)
            if (ctx.A[ctx.r_indx[ci]].size() == 0)
                ctx.Vi.diagonal().array() += prev_varcmp[ci];

        // LLT-only path (factorize_only, no diagV_adj)
        if (factorize_only && !ctx.reml_diagV_adj) {
            gcta_blas_int blas_n = static_cast<gcta_blas_int>(ctx.n);
            if (gcta_dpotrf(blas_n, ctx.Vi.data(), blas_n) == 0) {
                logdet = 2.0 * ctx.Vi.diagonal().array().log().sum();
                ctx.Vi_L.swap(ctx.Vi);
                ctx.Vi.resize(0, 0);
                ctx.Vi_use_llt = true;
                return true;
            }
            // dpotrf failed: reassemble
            ctx.Vi.triangularView<Eigen::Lower>().setZero();
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < ctx.n; j++)
                for (int ci = 0; ci < num_comp; ci++)
                    if (ctx.A[ctx.r_indx[ci]].size() > 0)
                        ctx.Vi.col(j).tail(ctx.n - j) +=
                            prev_varcmp[ci] * ctx.A[ctx.r_indx[ci]].col(j).tail(ctx.n - j);
            for (int ci = 0; ci < num_comp; ci++)
                if (ctx.A[ctx.r_indx[ci]].size() == 0)
                    ctx.Vi.diagonal().array() += prev_varcmp[ci];
        }

        INVmethod method_try = ctx.reml_inv_mtd ? static_cast<INVmethod>(ctx.reml_inv_mtd) : INV_LLT;
        int rank = 0;
        bool ret = true;

        if (method_try == INV_LLT && !factorize_only) {
            gcta_blas_int blas_n_f = static_cast<gcta_blas_int>(ctx.n);
            bool llt_ok = (gcta_dpotrf(blas_n_f, ctx.Vi.data(), blas_n_f) == 0);
            if (llt_ok) {
                logdet = ctx.Vi.diagonal().array().square().log().sum();
                llt_ok = (gcta_dpotri(blas_n_f, ctx.Vi.data(), blas_n_f) == 0);
            }
            if (llt_ok) {
                ctx.Vi.triangularView<Eigen::Upper>() = ctx.Vi.transpose();
                return true;
            }
            // Reassemble for LU fallback
            ctx.Vi.triangularView<Eigen::Lower>().setZero();
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < ctx.n; j++)
                for (int ci = 0; ci < num_comp; ci++)
                    if (ctx.A[ctx.r_indx[ci]].size() > 0)
                        ctx.Vi.col(j).tail(ctx.n - j) +=
                            prev_varcmp[ci] * ctx.A[ctx.r_indx[ci]].col(j).tail(ctx.n - j);
            for (int ci = 0; ci < num_comp; ci++)
                if (ctx.A[ctx.r_indx[ci]].size() == 0)
                    ctx.Vi.diagonal().array() += prev_varcmp[ci];
            ctx.Vi.triangularView<Eigen::Upper>() =
                ctx.Vi.triangularView<Eigen::Lower>().transpose();
            method_try = INV_LU;
        }

        if (!SquareMatrixInverse(ctx.Vi, logdet, rank, method_try)) {
            LOGGER << "Warning: V matrix is not invertible." << std::endl;
            if (ctx.reml_diagV_adj == 1) {
                LOGGER << "A small positive value is added to the diagonals. Results might not be reliable!" << std::endl;
                double d_buf = ctx.Vi.diagonal().mean() * ctx.reml_diag_mul;
                for (int j = 0; j < ctx.n; j++) ctx.Vi(j, j) += d_buf;
                if (!SquareMatrixInverse(ctx.Vi, logdet, rank, method_try)) {
                    LOGGER << "Still can't be inverted." << std::endl;
                    ret = false;
                }
            } else if (ctx.reml_diagV_adj == 2) {
                LOGGER << "Switching to the bending approach." << std::endl;
                bend_V(ctx, ctx.Vi);
            } else {
                LOGGER.e(0, "the variance-covariance matrix V is not invertible.");
                ret = false;
            }
        }
        return ret;
    }
    return true;
}

double calcu_P_impl(RemlCtx& ctx, RemlMat* P) {
    // Build Vi_X = V^{-1} X
    if (ctx.Vi_use_woodbury) {
        if (ctx.UkTX.rows() != ctx.woodbury_rank_ || ctx.UkTX.cols() != ctx.X_c)
            ctx.UkTX.noalias() = ctx.Uk.transpose() * ctx.X;
        if ((int)ctx.UkTy.size() != ctx.woodbury_rank_)
            ctx.UkTy.noalias() = ctx.Uk.transpose() * ctx.y;
        const RemlMat ck_UkTX = ctx.ck.asDiagonal() * ctx.UkTX;
        ctx.Vi_X = (ctx.X - ctx.Uk * ck_UkTX) / ctx.sigma2_eff;
        ctx.Uk_Vi_X = (ctx.UkTX - ck_UkTX) / ctx.sigma2_eff;
    } else if (ctx.Vi_use_llt) {
        ctx.Vi_X = ctx.X;
        ctx.Vi_L.triangularView<Eigen::Lower>().solveInPlace(ctx.Vi_X);
        ctx.Vi_L.triangularView<Eigen::Lower>().adjoint().solveInPlace(ctx.Vi_X);
    } else {
        ctx.Vi_X.noalias() = ctx.Vi.selfadjointView<Eigen::Lower>() * ctx.X;
    }
    ctx.Xt_Vi_X_i.noalias() = ctx.X.transpose() * ctx.Vi_X;

    double logdet_Xt_Vi_X = 0.0;
    int rank = 0;
    INVmethod method = (ctx.reml_inv_mtd == 0) ? INV_LLT : static_cast<INVmethod>(ctx.reml_inv_mtd);
    if (!SquareMatrixInverse(ctx.Xt_Vi_X_i, logdet_Xt_Vi_X, rank, method))
        LOGGER.e(0, "the X'V^{-1}X matrix is not invertible. Please check covariates.");

    if (!P) return logdet_Xt_Vi_X;

    // Materialise V^{-1} if we need the full P matrix
    if (ctx.Vi_use_woodbury) {
        ctx.Vi.resize(ctx.n, ctx.n);
        ctx.Vi.setIdentity();
        ctx.Vi /= ctx.sigma2_eff;
        RemlMat Uk_scaled = ctx.Uk;
        for (int j = 0; j < ctx.woodbury_rank_; ++j)
            Uk_scaled.col(j) *= std::sqrt(ctx.ck[j] / ctx.sigma2_eff);
        ctx.Vi.selfadjointView<Eigen::Lower>().rankUpdate(Uk_scaled, -1.0);
        ctx.Vi.triangularView<Eigen::Upper>() = ctx.Vi.transpose();
    } else if (ctx.Vi_use_llt) {
        ctx.Vi.swap(ctx.Vi_L);
        gcta_blas_int blas_n_p = static_cast<gcta_blas_int>(ctx.n);
        if (gcta_dpotri(blas_n_p, ctx.Vi.data(), blas_n_p) != 0)
            LOGGER.e(0, "dpotri failed when materialising V^{-1} for P.");
        ctx.Vi.triangularView<Eigen::Upper>() = ctx.Vi.transpose();
        ctx.Vi_use_llt = false;
    }

    Eigen::LLT<RemlMat> llt(ctx.Xt_Vi_X_i);
    if (llt.info() == Eigen::Success) {
        RemlMat Z;
        Z.noalias() = ctx.Vi_X * llt.matrixL();
        P->swap(ctx.Vi);
        P->selfadjointView<Eigen::Lower>().rankUpdate(Z, -1.0);
        P->triangularView<Eigen::Upper>() = P->transpose();
    } else {
        RemlMat W;
        W.noalias() = ctx.Vi_X * ctx.Xt_Vi_X_i;
        P->swap(ctx.Vi);
        P->noalias() -= W * ctx.Vi_X.transpose();
    }
    return logdet_Xt_Vi_X;
}

void calcu_tr_PA_woodbury(const RemlCtx& ctx, RemlVec& tr_PA) {
    const int ncomp = static_cast<int>(ctx.r_indx.size());
    tr_PA.resize(ncomp);

    const double sigma2_eff = ctx.sigma2_eff;
    const double sg2        = ctx.sg2;
    const double se2        = sigma2_eff - sg2 * ctx.lambda_tail;
    const double lambda_t   = ctx.lambda_tail;
    const int    n          = ctx.n;
    const int    k          = ctx.woodbury_rank_;
    const double sum_ck     = ctx.ck.sum();
    const double tr_Vinv    = (n - sum_ck) / sigma2_eff;
    const double tr_Vinv_K  = (sg2 > 1e-15)
        ? (n * lambda_t / sigma2_eff + (se2 / (sg2 * sigma2_eff)) * sum_ck)
        : tr_Vinv * lambda_t;

    const int c = static_cast<int>(ctx.Vi_X.cols());
    RemlMat ViXTViX(c, c);
    ViXTViX.noalias() = ctx.Vi_X.transpose() * ctx.Vi_X;
    const double tr_corr_I = (ctx.Xt_Vi_X_i * ViXTViX).trace();

    tr_PA(ncomp - 1) = tr_Vinv - tr_corr_I;

    RemlVec delta = ctx.dk.array() - ctx.lambda_tail;
    RemlMat delta_UkViX = delta.asDiagonal() * ctx.Uk_Vi_X;
    RemlMat ViX_K_ViX = lambda_t * ViXTViX;
    ViX_K_ViX.noalias() += ctx.Uk_Vi_X.transpose() * delta_UkViX;
    const double tr_corr_K = (ctx.Xt_Vi_X_i * ViX_K_ViX).trace();

    tr_PA(0) = tr_Vinv_K - tr_corr_K;
}

void calcu_tr_PA_hutchpp(RemlCtx& ctx, RemlVec& tr_PA, int m_probes) {
    const int ncomp = static_cast<int>(ctx.r_indx.size());
    tr_PA.resize(ncomp);
    const int k = std::max(m_probes / 3, 3);

    if (ctx.hutchpp_S.rows() != ctx.n || ctx.hutchpp_S.cols() != k) {
        const uint32_t seed = static_cast<uint32_t>(ctx.n) * 2654435761u
                              ^ (static_cast<uint32_t>(ncomp) * 2246822519u)
                              ^ 0x9e3779b9u;
        std::mt19937 probe_rng(seed);
        std::uniform_int_distribution<int> coin(0, 1);
        ctx.hutchpp_S.resize(ctx.n, k);
        ctx.hutchpp_G.resize(ctx.n, k);
        for (int j = 0; j < k; j++)
            for (int r = 0; r < ctx.n; r++) {
                ctx.hutchpp_S(r, j) = coin(probe_rng) ? 1.0 : -1.0;
                ctx.hutchpp_G(r, j) = coin(probe_rng) ? 1.0 : -1.0;
            }
    }

    for (int ci = 0; ci < ncomp; ci++) {
        const bool is_I = (ctx.A[ctx.r_indx[ci]].size() == 0);

        auto applyPA_mat = [&](const RemlMat& Z) -> RemlMat {
            if (ctx.Vi_use_woodbury && ci == 0) {
                RemlMat UkZ = ctx.Uk.transpose() * Z;
                UkZ.array().colwise() *= (ctx.dk.array() - ctx.lambda_tail);
                RemlMat KZ = ctx.lambda_tail * Z + ctx.Uk * UkZ;
                return applyP_mat(ctx, KZ);
            }
            return is_I ? applyP_mat(ctx, Z)
                        : applyP_mat(ctx, RemlMat(ctx.A[ctx.r_indx[ci]].selfadjointView<Eigen::Lower>() * Z));
        };

        RemlMat K = applyPA_mat(ctx.hutchpp_S);
        for (int pw = 0; pw < ctx.reml_trace_power_iter; pw++) {
            Eigen::HouseholderQR<RemlMat> qr_pw(K);
            K = qr_pw.householderQ() * RemlMat::Identity(ctx.n, k);
            K = applyPA_mat(K);
        }
        Eigen::HouseholderQR<RemlMat> qr(K);
        RemlMat Q = qr.householderQ() * RemlMat::Identity(ctx.n, k);

        const RemlMat MQ = applyPA_mat(Q);
        const double t_lr = Q.cwiseProduct(MQ).sum();
        const RemlMat MG = applyPA_mat(ctx.hutchpp_G);

        const RemlMat QtG = Q.transpose() * ctx.hutchpp_G;
        const RemlMat R_  = ctx.hutchpp_G - Q * QtG;
        const RemlMat MR  = MG - MQ * QtG;

        tr_PA(ci) = t_lr + R_.cwiseProduct(MR).sum() / k;
    }
}

void calcu_tr_PA(const RemlCtx& ctx, const RemlMat& P, RemlVec& tr_PA) {
    const int m = static_cast<int>(ctx.r_indx.size());
    tr_PA.resize(m);
    for (int i = 0; i < m; i++) {
        if (ctx.A[ctx.r_indx[i]].size() == 0) {
            double s = 0.0;
            #pragma omp parallel for reduction(+:s) schedule(static)
            for (int col = 0; col < ctx.n; col++) s += P(col, col);
            tr_PA(i) = s;
        } else {
            const auto& Ai = ctx.A[ctx.r_indx[i]];
            double s = 0.0;
            #pragma omp parallel for reduction(+:s) schedule(guided)
            for (int col = 0; col < ctx.n; col++) {
                const int tail = ctx.n - col - 1;
                s += P(col, col) * Ai(col, col);
                if (tail > 0)
                    s += 2.0 * P.col(col).tail(tail).dot(Ai.col(col).tail(tail));
            }
            tr_PA(i) = s;
        }
    }
}

void calcu_Hi(RemlCtx& ctx, RemlMat& P, RemlMat& Hi) {
    P = (P + P.transpose()) * 0.5;
    const int m = static_cast<int>(ctx.r_indx.size());

    std::vector<RemlMat> PA(m);
    for (int i = 0; i < m; i++) {
        if (ctx.Vi_use_woodbury && i == 0) {
            RemlMat PUk = P.selfadjointView<Eigen::Upper>() * ctx.Uk;
            RemlVec delta_v = ctx.dk.array() - ctx.lambda_tail;
            PA[i] = ctx.lambda_tail * P;
            PA[i].noalias() += PUk * delta_v.asDiagonal() * ctx.Uk.transpose();
        } else if (ctx.A[ctx.r_indx[i]].size() == 0) {
            PA[i].resize(0, 0);
        } else {
            PA[i].noalias() = P.selfadjointView<Eigen::Upper>() * ctx.A[ctx.r_indx[i]];
        }
    }

    for (int i = 0; i < m; i++) {
        const bool i_id = (PA[i].size() == 0);
        for (int j = 0; j <= i; j++) {
            const bool j_id = (PA[j].size() == 0);
            double val;
            if (i_id && j_id)       val = P.squaredNorm();
            else if (i_id)          val = PA[j].cwiseProduct(P).sum();
            else if (j_id)          val = PA[i].cwiseProduct(P).sum();
            else                    val = PA[i].cwiseProduct(PA[j].transpose()).sum();
            Hi(i, j) = Hi(j, i) = val;
        }
    }

    if (!inverse_H(ctx, Hi)) {
        if (ctx.reml_force_converge) {
            LOGGER.w(0, "the information matrix is not invertible.");
            ctx.reml_AI_not_invertible = true;
        } else {
            LOGGER.e(0, "the information matrix is not invertible.");
        }
    }
}

void ai_reml(RemlCtx& ctx, RemlMat& P, RemlMat& Hi, RemlVec& Py,
             RemlVec& prev_varcmp, RemlVec& varcmp, double dlogL) {
    const bool use_approx     = ctx.reml_trace_approx;
    const bool woodbury_active = ctx.Vi_use_woodbury;

    if (use_approx || woodbury_active)
        Py = applyP_vec(ctx, ctx.y);
    else
        Py.noalias() = P.selfadjointView<Eigen::Lower>() * ctx.y;

    const int m = static_cast<int>(ctx.r_indx.size());
    RemlMat APy(ctx.n, m);
    for (int i = 0; i < m; i++) {
        if (woodbury_active && i == 0)
            APy.col(i) = woodbury_Kv(ctx, Py);
        else if (ctx.A[ctx.r_indx[i]].size() == 0)
            APy.col(i) = Py;
        else
            APy.col(i).noalias() = ctx.A[ctx.r_indx[i]].selfadjointView<Eigen::Lower>() * Py;
    }

    RemlVec R(m);
    if (use_approx || woodbury_active) {
        R.noalias() = APy.transpose() * Py;
        const RemlMat PAPy = applyP_mat(ctx, APy);
        Hi.noalias() = APy.transpose() * PAPy;
    } else {
        R.noalias() = APy.transpose() * Py;
        RemlMat PAPy(ctx.n, m);
        PAPy.noalias() = P.selfadjointView<Eigen::Lower>() * APy;
        Hi.noalias() = APy.transpose() * PAPy;
    }
    Hi = 0.5 * Hi;

    RemlVec tr_PA;
    if (woodbury_active) {
        calcu_tr_PA_woodbury(ctx, tr_PA);
    } else if (use_approx) {
        const int eff_nprobes = (std::fabs(dlogL) < 1.0)
            ? ctx.reml_trace_approx_nprobes
            : std::max(ctx.reml_trace_approx_nprobes / 3, 9);
        calcu_tr_PA_hutchpp(ctx, tr_PA, eff_nprobes);
    } else {
        calcu_tr_PA(ctx, P, tr_PA);
    }
    R = -0.5 * (tr_PA - R);

    if (!inverse_H(ctx, Hi)) {
        if (ctx.reml_force_converge) {
            LOGGER.w(0, "the information matrix is not invertible.");
            ctx.reml_AI_not_invertible = true;
            return;
        } else {
            LOGGER.e(0, "the information matrix is not invertible.");
        }
    }
    RemlVec delta = Hi * R;
    if (dlogL > 1.0) varcmp = prev_varcmp + 0.316 * delta;
    else             varcmp = prev_varcmp + delta;
}

void em_reml(RemlCtx& ctx, RemlMat& P, RemlVec& Py,
             RemlVec& prev_varcmp, RemlVec& varcmp, double dlogL) {
    const bool woodbury_active = ctx.Vi_use_woodbury;
    const bool use_approx     = ctx.reml_trace_approx;

    RemlVec tr_PA;
    if (woodbury_active) {
        calcu_tr_PA_woodbury(ctx, tr_PA);
        Py = applyP_vec(ctx, ctx.y);
    } else if (use_approx) {
        const int eff_nprobes = (std::fabs(dlogL) < 1.0)
            ? ctx.reml_trace_approx_nprobes
            : std::max(ctx.reml_trace_approx_nprobes / 3, 9);
        calcu_tr_PA_hutchpp(ctx, tr_PA, eff_nprobes);
        Py = applyP_vec(ctx, ctx.y);
    } else {
        calcu_tr_PA(ctx, P, tr_PA);
        Py.noalias() = P.selfadjointView<Eigen::Lower>() * ctx.y;
    }

    const int m = static_cast<int>(ctx.r_indx.size());
    RemlVec R(m);
    RemlVec tmp(ctx.n);
    for (int i = 0; i < m; i++) {
        if (woodbury_active && i == 0) {
            tmp = woodbury_Kv(ctx, Py);
            R(i) = Py.dot(tmp);
        } else if (ctx.A[ctx.r_indx[i]].size() == 0) {
            R(i) = Py.squaredNorm();
        } else {
            tmp.noalias() = ctx.A[ctx.r_indx[i]].selfadjointView<Eigen::Lower>() * Py;
            R(i) = Py.dot(tmp);
        }
        varcmp(i) = prev_varcmp(i) - prev_varcmp(i) * prev_varcmp(i) * (tr_PA(i) - R(i)) / ctx.n;
    }
}

double reml_iteration(RemlCtx& ctx,
                      RemlMat& Vi_X_out, RemlMat& Xt_Vi_X_i_out, RemlMat& Hi,
                      RemlVec& Py, RemlVec& varcmp,
                      bool prior_var_flag, bool no_constrain) {
    std::vector<std::string> mtd_str = {"AI-REML", "Fisher-scoring REML", "EM-REML"};
    const int m = static_cast<int>(ctx.r_indx.size());
    int constrain_num = 0, iter = 0;
    int reml_mtd_tmp = ctx.reml_mtd;
    double logdet = 0.0, logdet_Xt_Vi_X = 0.0;
    double prev_lgL = -1e20, lgL = -1e20, dlogL = 1000.0;
    RemlVec prev_prev_varcmp(varcmp), prev_varcmp(varcmp), varcomp_init(varcmp);

    if (ctx.reml_trace_approx && !ctx.Vi_use_woodbury) {
        LOGGER << "Using Hutch++ stochastic trace estimator with "
               << ctx.reml_trace_approx_nprobes << " probes." << std::endl;
    }

    bool converged_flag = false;
    for (iter = 0; iter < ctx.reml_max_iter; iter++) {
        if (iter == 0) {
            prev_varcmp = varcomp_init;
            if (prior_var_flag) {
                if (ctx.reml_fixed_var)
                    LOGGER << "Variance components fixed at: " << varcmp.transpose() << std::endl;
                else
                    LOGGER << "Prior values of variance components: " << varcmp.transpose() << std::endl;
            } else {
                ctx.reml_mtd = 2;
                LOGGER << "Calculating prior values of variance components by EM-REML ..." << std::endl;
            }
        }
        if (iter == 1) {
            ctx.reml_mtd = reml_mtd_tmp;
            LOGGER << "Running " << mtd_str[ctx.reml_mtd] << " algorithm ...\nIter.\tlogL\t";
            for (int i = 0; i < m; i++) LOGGER << ctx.var_name[ctx.r_indx[i]] << "\t";
            LOGGER << std::endl;
        }

        const bool skip_P_this_iter = (ctx.reml_trace_approx || ctx.Vi_use_woodbury) && ctx.reml_mtd != 1;
        if (!calcu_Vi(ctx, prev_varcmp, logdet, iter, skip_P_this_iter)) {
            LOGGER << "Warning: V matrix is not positive-definite." << std::endl;
            varcmp = prev_prev_varcmp;
            if (!calcu_Vi(ctx, varcmp, logdet, iter, false))
                LOGGER.e(0, "V matrix is not positive-definite.");
            // Rebuild P for Hi computation in break path
            RemlMat P_tmp(ctx.n, ctx.n);
            calcu_P_impl(ctx, &P_tmp);
            calcu_Hi(ctx, P_tmp, Hi);
            Hi = 2 * Hi;
            break;
        }

        double logdet_Xt_Vi_X2;
        if (skip_P_this_iter) {
            logdet_Xt_Vi_X2 = calcu_P_impl(ctx, nullptr);
            ctx.P.resize(0, 0);
        } else {
            logdet_Xt_Vi_X2 = calcu_P_impl(ctx, &ctx.P);
            ctx.Vi.resize(0, 0);
        }
        logdet_Xt_Vi_X = logdet_Xt_Vi_X2;

        if (ctx.reml_mtd == 0)
            ai_reml(ctx, ctx.P, Hi, Py, prev_varcmp, varcmp, dlogL);
        else if (ctx.reml_mtd == 2)
            em_reml(ctx, ctx.P, Py, prev_varcmp, varcmp, dlogL);
        // mtd 1 (Fisher scoring) requires full P — not supported in skip_P path
        // but skip_P_this_iter is false for mtd 1, so this is safe.

        lgL = -0.5 * (logdet_Xt_Vi_X + logdet + (ctx.y.transpose() * Py)(0, 0));

        if (ctx.reml_force_converge && ctx.reml_AI_not_invertible) break;

        if (!no_constrain) constrain_num = constrain_varcmp(ctx, varcmp);

        if (iter > 0) {
            LOGGER << iter << "\t" << std::fixed << LOGGER.setprecision(2) << lgL << "\t";
            for (int i = 0; i < m; i++) LOGGER << LOGGER.setprecision(5) << varcmp[i] << "\t";
            if (constrain_num > 0) LOGGER << "(" << constrain_num << " component(s) constrained)" << std::endl;
            else LOGGER << std::endl;
        } else {
            if (!prior_var_flag) LOGGER << "Updated prior values: " << varcmp.transpose() << std::endl;
            LOGGER << "logL: " << lgL << std::endl;
        }

        if (ctx.reml_fixed_var) { varcmp = prev_varcmp; break; }

        if (constrain_num * 2 > m) {
            if (ctx.reml_allow_constrain_run) {
                LOGGER.w(0, "more than half of the variance components are constrained.");
            } else {
                LOGGER.e(0, "analysis stopped because more than half of the variance components are constrained.");
            }
        }

        if ((ctx.reml_force_converge || ctx.reml_no_converge) && prev_lgL > lgL) {
            varcmp = prev_varcmp;
            RemlMat P_tmp;
            calcu_P_impl(ctx, &P_tmp);
            calcu_Hi(ctx, P_tmp, Hi);
            Hi = 2 * Hi;
            break;
        }

        dlogL = lgL - prev_lgL;
        if ((varcmp - prev_varcmp).squaredNorm() / varcmp.squaredNorm() < 1e-8
                && (std::fabs(dlogL) < 1e-4 || (std::fabs(dlogL) < 1e-2 && dlogL < 0))) {
            converged_flag = true;
            if (ctx.reml_mtd == 2) {
                RemlMat P_tmp;
                calcu_P_impl(ctx, &P_tmp);
                calcu_Hi(ctx, P_tmp, Hi);
                Hi = 2 * Hi;
            }
            break;
        }

        ctx.Vi.swap(ctx.P);
        prev_prev_varcmp = prev_varcmp;
        prev_varcmp = varcmp;
        prev_lgL = lgL;
    }

    if (ctx.reml_fixed_var)
        LOGGER.w(0, "model evaluated at fixed variance components; log-likelihood may not be maximised.");
    else {
        if (converged_flag) LOGGER << "Log-likelihood ratio converged." << std::endl;
        else if (ctx.reml_force_converge || ctx.reml_no_converge)
            LOGGER.w(0, "Log-likelihood not converged. Results are not reliable.");
        else if (iter == ctx.reml_max_iter) {
            std::ostringstream errmsg;
            errmsg << "Log-likelihood not converged (stop after " << ctx.reml_max_iter << " iterations). "
                   << "Use --reml-maxit for more iterations.";
            if (ctx.reml_max_iter > 1) LOGGER.e(0, errmsg.str());
        }
    }
    ctx.P.resize(0, 0);

    // Export final Vi_X and Xt_Vi_X_i to output parameters
    Vi_X_out    = ctx.Vi_X;
    Xt_Vi_X_i_out = ctx.Xt_Vi_X_i;

    return lgL;
}

} // anonymous namespace

// ─────────────────────────────────────────────────────────────────────────────
// Public API implementations
// ─────────────────────────────────────────────────────────────────────────────

namespace reml {

void compute_woodbury_basis(RemlCtx& ctx) {
    if (ctx.reml_mtd == 1)
        LOGGER.e(0, "--reml-woodbury is incompatible with Fisher-scoring REML.");
    if ((int)ctx.r_indx.size() != 2)
        LOGGER.e(0, "--reml-woodbury supports only single-GRM models.");
    if (ctx.A[ctx.r_indx[0]].size() == 0)
        LOGGER.e(0, "--reml-woodbury: GRM component is identity; cannot compute basis.");

    const bool auto_k = (ctx.woodbury_buffer_factor > 0.0);
    const int  n      = ctx.n;
    int k = ctx.woodbury_rank;  // fixed-k or starting point; overridden in auto-k

    int k_svd;
    if (auto_k) {
        k_svd = (ctx.woodbury_k_max > 0) ? ctx.woodbury_k_max : std::min(n - 1, 1200);
        LOGGER << "\nComputing Woodbury basis (auto-k, k_max=" << k_svd
               << ", buffer=" << ctx.woodbury_buffer_factor << ") ..." << std::endl;
    } else {
        k_svd = k;
        LOGGER << "\nComputing Woodbury low-rank basis (k=" << k << ") ..." << std::endl;
    }
    if (k_svd >= n) LOGGER.e(0, "--reml-woodbury rank must be < n.");

    const Eigen::MatrixXd& K_dbl = ctx.A[ctx.r_indx[0]];

    const int oversample = 20;
    const int k_ext = std::min(k_svd + oversample, n - 1);
    Eigen::MatrixXd omega;
    {
        const int k_prev = static_cast<int>(ctx.Uk.cols());
        const bool has_warm = (!ctx.woodbury_nystrom && k_prev > 0 && ctx.Uk.rows() == n);
        if (has_warm) {
            omega.resize(n, k_ext);
            const int k_copy = std::min(k_prev, k_ext);
            omega.leftCols(k_copy)  = ctx.Uk.leftCols(k_copy);
            if (k_copy < k_ext)
                omega.rightCols(k_ext - k_copy) = Eigen::MatrixXd::Random(n, k_ext - k_copy);
        } else {
            omega = Eigen::MatrixXd::Random(n, k_ext);
        }
    }
    Eigen::MatrixXd Y = K_dbl.selfadjointView<Eigen::Lower>() * omega;

    Eigen::VectorXd eval_full;
    Eigen::MatrixXd evec_full;

    if (ctx.woodbury_nystrom) {
        Eigen::MatrixXd C = omega.transpose() * Y;
        omega.resize(0, 0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_C(C);
        if (es_C.info() != Eigen::Success)
            LOGGER.e(0, "Woodbury Nystrom: eigendecomposition of sketch C failed.");
        const double lam_max = es_C.eigenvalues().maxCoeff();
        const double eps_C   = 1e-8 * std::max(lam_max, 1.0);
        Eigen::VectorXd lam_sqrt_inv =
            es_C.eigenvalues().cwiseMax(eps_C).cwiseSqrt().cwiseInverse();
        Y = Y * (es_C.eigenvectors() * lam_sqrt_inv.asDiagonal());
        Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU> svd(Y);
        Y.resize(0, 0);
        eval_full = svd.singularValues().head(k_svd).array().square();
        evec_full = svd.matrixU().leftCols(k_svd);
    } else {
        omega.resize(0, 0);
        constexpr int power_iter = 3;
        for (int pi = 0; pi < power_iter; ++pi) {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr_pi(Y);
            Y = qr_pi.householderQ() * Eigen::MatrixXd::Identity(n, k_ext);
            Y = K_dbl.selfadjointView<Eigen::Lower>() * Y;
        }
        {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(Y);
            Y = qr.householderQ() * Eigen::MatrixXd::Identity(n, k_ext);
        }
        Eigen::MatrixXd KY = K_dbl.selfadjointView<Eigen::Lower>() * Y;
        Eigen::MatrixXd B  = Y.transpose() * KY;
        KY.resize(0, 0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
        if (es.info() != Eigen::Success)
            LOGGER.e(0, "Woodbury: eigendecomposition of projected GRM failed.");
        eval_full = es.eigenvalues().tail(k_svd).reverse();
        evec_full = Y * es.eigenvectors().rightCols(k_svd).rowwise().reverse();
    }

    if (auto_k) {
        double M = 0.0;
        if (ctx.grm_N.rows() == n && ctx.grm_N.cols() == n)
            M = ctx.grm_N.diagonal().mean();
        if (M <= 0.0)
            LOGGER.e(0, "--reml-woodbury auto: cannot determine SNP count. Use --reml-woodbury <k>.");
        const double gamma       = static_cast<double>(n) / M;
        const double lambda_plus = std::pow(1.0 + std::sqrt(gamma), 2.0);
        const int    k_signal   = static_cast<int>((eval_full.array() > lambda_plus).count());
        const int    k_buffered = std::max(20, static_cast<int>(std::ceil(k_signal * ctx.woodbury_buffer_factor)));
        k = std::min(k_svd, k_buffered);
        LOGGER << "MP bulk edge lambda+ = " << lambda_plus
               << " (n=" << n << ", M=" << static_cast<long long>(M) << ")"
               << ", eigenvalues above lambda+ = " << k_signal
               << ", using k = " << k << " (buffer=" << ctx.woodbury_buffer_factor << ")" << std::endl;
        if (k_buffered > k_svd)
            LOGGER.w(0, "Woodbury auto-k: buffer-implied k=" + std::to_string(k_buffered)
                     + " exceeds k_max=" + std::to_string(k_svd) + "; clamped to k_max.");
    }

    Eigen::VectorXd eval = eval_full.head(k);
    Eigen::MatrixXd evec = evec_full.leftCols(k);

    const double trace_K = K_dbl.diagonal().sum();
    ctx.lambda_tail = (trace_K - eval.sum()) / static_cast<double>(n - k);
    if (ctx.lambda_tail < 0.0) {
        LOGGER.w(0, "Woodbury: lambda_tail < 0 (" + std::to_string(ctx.lambda_tail) + "); clamped to 0.");
        ctx.lambda_tail = 0.0;
    }

    {
        double diag_sq = K_dbl.diagonal().squaredNorm();
        double off_sq  = 0.0;
        #pragma omp parallel for reduction(+:off_sq) schedule(static)
        for (int j = 0; j < n; ++j)
            off_sq += K_dbl.col(j).tail(n - j - 1).squaredNorm();
        const double trace_K2    = diag_sq + 2.0 * off_sq;
        const double tail_sum_sq = trace_K2 - eval.squaredNorm();
        ctx.tail_d_var = std::max(0.0, tail_sum_sq / static_cast<double>(n - k)
                                        - ctx.lambda_tail * ctx.lambda_tail);
    }

    eval = eval.cwiseMax(0.0);
    ctx.dk            = eval;
    ctx.Uk            = evec;
    ctx.woodbury_rank_ = k;

    ctx.A[ctx.r_indx[0]].resize(0, 0);
    ctx.Vi_use_woodbury = true;

    LOGGER << "Woodbury basis: k=" << k
           << ", lambda_tail=" << ctx.lambda_tail
           << ", tail_d_var=" << ctx.tail_d_var << std::endl;
}

void compute(RemlCtx& ctx,
             const std::vector<double>& priors,
             const std::vector<double>& priors_var,
             bool no_constrain) {
    const bool priors_flag = !priors.empty() || !priors_var.empty();

    {
        RemlVec y_tmp = ctx.y.array() - ctx.y.mean();
        ctx.y_Ssq = y_tmp.squaredNorm() / (ctx.n - 1.0);
        if (!(std::fabs(ctx.y_Ssq) < 1e30))
            LOGGER.e(0, "phenotypic variance is infinite. Check phenotype file.");
    }

    if (priors_flag && (int)priors_var.size() < (int)ctx.r_indx.size() - 1) {
        std::ostringstream errmsg;
        errmsg << "in option --reml-priors-var. There are " << ctx.r_indx.size()
               << " variance components. At least " << ctx.r_indx.size() - 1
               << " prior values should be specified.";
        LOGGER.e(0, errmsg.str());
    }

    LOGGER << "\nPerforming REML analysis ... (Note: may take hours depending on sample size)." << std::endl;
    if (ctx.n < 10) LOGGER.e(0, "sample size is too small.");
    LOGGER << ctx.n << " observations, " << ctx.X_c << " fixed effect(s), and "
           << ctx.r_indx.size() << " variance component(s) (including residual)." << std::endl;

    // Woodbury basis (must be done before init_varcomp for HE warm-start)
    if (ctx.woodbury_rank > 0)
        compute_woodbury_basis(ctx);
    else if (ctx.woodbury_rank < 0) {
        // auto-k: temporarily set rank to 0, let compute_woodbury_basis use buffer_factor
        ctx.woodbury_rank = 0; // auto-k triggered by buffer_factor > 0
        compute_woodbury_basis(ctx);
    }

    RemlMat Vi_X_out(ctx.n, ctx.X_c), Xt_Vi_X_i_out(ctx.X_c, ctx.X_c);
    RemlMat Hi(ctx.r_indx.size(), ctx.r_indx.size());
    RemlVec Py(ctx.n);
    RemlVec varcmp;
    init_varcomp(ctx, priors_var, priors, varcmp);

    /*double lgL =*/ reml_iteration(ctx, Vi_X_out, Xt_Vi_X_i_out, Hi, Py, varcmp, priors_flag, no_constrain);

    // Compute fixed effects: b = (X'V^{-1}X)^{-1} X'V^{-1}y
    ctx.b = Xt_Vi_X_i_out * (Vi_X_out.transpose() * ctx.y);

    // Store variance components as a std::vector
    ctx.varcmp.resize(ctx.r_indx.size());
    for (int i = 0; i < (int)ctx.r_indx.size(); i++) ctx.varcmp[i] = varcmp[i];

    // Ensure V^{-1} is available for downstream MLMA streaming.
    // If neither Woodbury nor Hutch++ path: calcu_Vi with factorize_only=true so
    // we get L in ctx.Vi_L (or Vi in ctx.Vi) without needing P.
    if (!ctx.reml_trace_approx && !ctx.Vi_use_woodbury) {
        double logdet_dummy = 0.0;
        int iter_dummy = 0;
        calcu_Vi(ctx, varcmp, logdet_dummy, iter_dummy, /*factorize_only=*/true);
    }
}

RemlState build_reml_state(const RemlCtx& ctx) {
    RemlState rs;
    rs.n           = static_cast<int32_t>(ctx.n);
    rs.x_c         = static_cast<int32_t>(ctx.X_c);
    rs.is_woodbury = ctx.Vi_use_woodbury;

    // Fixed-effects vector (double → float)
    rs.b = ctx.b.cast<float>();

    if (ctx.Vi_use_woodbury) {
        // WoodburyMLMACache stores Uk_f as k×n (transposed) for GEMM efficiency.
        // ctx.Uk is n×k — transpose before storing.
        rs.wb.Uk_f         = ctx.Uk.transpose().cast<float>();  // k×n
        rs.wb.ck_f         = ctx.ck.cast<float>();
        rs.wb.sqrt_ck_f    = rs.wb.ck_f.cwiseSqrt();
        rs.wb.sigma2_eff_f = static_cast<float>(ctx.sigma2_eff);
    } else if (ctx.Vi_use_llt) {
        // Vi_L holds the Cholesky factor L; materialise V^{-1} via dpotri.
        // Copy to avoid mutating ctx (it's const ref).
        Eigen::MatrixXd Vi_copy = ctx.Vi_L;
        gcta_blas_int blas_n = static_cast<gcta_blas_int>(ctx.n);
        if (gcta_dpotri(blas_n, Vi_copy.data(), blas_n) != 0)
            LOGGER.e(0, "build_reml_state: dpotri failed materialising V^{-1}.");
        Vi_copy.triangularView<Eigen::Upper>() = Vi_copy.transpose();
        rs.Vi = Vi_copy.cast<float>();
    } else {
        // Vi is already a full V^{-1}
        rs.Vi = ctx.Vi.cast<float>();
    }
    return rs;
}

} // namespace reml
