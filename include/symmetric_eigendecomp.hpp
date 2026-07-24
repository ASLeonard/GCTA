/*
 * symmetric_eigendecomp.hpp
 *
 * Approximate top-k eigendecomposition methods for a symmetric (typically
 * PSD) matrix A, shared between grm.cpp (gcta::pca) and RemlEngine.cpp
 * (compute_woodbury_basis).
 *
 * Every entry point takes a matvec functor rather than a concrete matrix
 * type: `apply(X)` must return A * X for X either an n-vector or an n x m
 * block. This lets callers pass a plain MatrixXd, a selfadjointView<Upper>
 * or selfadjointView<Lower> (grm.cpp's GRM is fully populated; RemlEngine's
 * ctx.A GRM only has the lower triangle valid), or in future a tiled/
 * streaming matvec, without this header caring.
 *
 * Methods:
 *   - randomized_symmetric_eigh  : Halko/Martinsson/Tropp randomized range
 *                                  finder + power iteration + Rayleigh-Ritz
 *                                  projection. Good default; k+oversample
 *                                  dense matvecs dominate cost.
 *   - lanczos_symmetric_eigh     : Spectra Lanczos (implicitly restarted
 *                                  Arnoldi on a symmetric operator), useful
 *                                  when the spectrum is well separated near
 *                                  the top and few matvecs are wanted.
 *   - tall_skinny_thin_svd       : QR-then-small-SVD reduction, used to
 *                                  extract the thin SVD of an n x k sketch
 *                                  (e.g. the Nystrom sketch) without paying
 *                                  for a bidiagonal divide-and-conquer SVD
 *                                  at n scale.
 */
#pragma once

#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace gcta_eigh {

struct EighResult {
    Eigen::VectorXd eigenvalues;   // descending, size k_target
    Eigen::MatrixXd eigenvectors;  // n x k_target, columns match eigenvalues
};

// Oversampling doesn't need to grow with k_target for generic top-k subspace
// accuracy (HMT's error bound doesn't require it). But callers that use
// eigenvalues near the *tail* of a large k_target block to locate a spectral
// threshold (e.g. Woodbury's auto-k Marchenko-Pastur edge, or its EIGMASS mass
// target) are relying on exactly the estimates that are least accurate —
// Rayleigh-Ritz quality degrades toward the edge of the requested block, and
// GRM eigenvalues cluster tightly there by construction (that's the point of
// the MP bulk edge). A fixed p=20 is a 3x buffer at k=6 (PCA; top eigenvalues
// are well-separated population-structure signal, not the bottleneck) but a
// 0.4% buffer at k=5000. Scale with k_target instead; the extra sketch
// columns are cheap once k_ext is already in the thousands.
inline int recommended_oversample(int k_target, int floor = 20, int cap = 200) {
    return std::clamp(k_target / 10, floor, cap);
}

// ─────────────────────────────────────────────────────────────────────────
// Randomized range finder + power iteration + Rayleigh-Ritz (rSVD)
// ─────────────────────────────────────────────────────────────────────────

// Build the initial randomized sketch Y = A * omega (n x k_ext).
// If `warm_start` is supplied (e.g. a previous Uk basis, or a PCA .eigenvec
// basis for the same GRM), its leading columns seed omega instead of a
// fresh Gaussian draw, which typically lets power_iterate_and_project()
// converge in fewer iterations.
template <typename MatVecApply>
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> build_randomized_sketch(
    MatVecApply&& apply,
    int n,
    int k_ext,
    const Eigen::MatrixXd* warm_start = nullptr)
{
    Eigen::MatrixXd omega;
    if (warm_start && warm_start->rows() == n && warm_start->cols() > 0) {
        omega.resize(n, k_ext);
        const int k_copy = std::min(static_cast<int>(warm_start->cols()), k_ext);
        omega.leftCols(k_copy) = warm_start->leftCols(k_copy);
        if (k_copy < k_ext)
            omega.rightCols(k_ext - k_copy) = Eigen::MatrixXd::Random(n, k_ext - k_copy);
    } else {
        omega = Eigen::MatrixXd::Random(n, k_ext);
    }
    Eigen::MatrixXd Y = apply(omega);
    return {std::move(omega), std::move(Y)};
}

// Core duplicated logic: given an initial sketch Y = A * omega (n x k_ext),
// re-orthogonalize/re-multiply for `power_iter` passes, then do a final QR
// and a Rayleigh-Ritz projection onto the k_ext x k_ext subspace to recover
// the k_target dominant eigenpairs of A.
//
// `apply` is called power_iter + 1 more times. Y is consumed/overwritten.
template <typename MatVecApply>
EighResult power_iterate_and_project(
    MatVecApply&& apply,
    Eigen::MatrixXd Y,
    int k_target,
    int power_iter = 3)
{
    const int n     = static_cast<int>(Y.rows());
    const int k_ext = static_cast<int>(Y.cols());
    if (k_target > k_ext)
        throw std::invalid_argument("power_iterate_and_project: k_target exceeds sketch width k_ext.");

    // Pre-allocated thin-Q scratch, reused across every QR in this call so
    // each pass avoids a fresh n x k_ext zero-init allocation (this matters:
    // at n=500k, k_ext~1200, each Identity(n,k_ext) is ~4.9 GB).
    Eigen::MatrixXd qr_scratch = Eigen::MatrixXd::Identity(n, k_ext);

    for (int pi = 0; pi < power_iter; ++pi) {
        Eigen::HouseholderQR<Eigen::MatrixXd> qr(Y);
        qr_scratch.setIdentity();
        Y = apply(qr.householderQ() * qr_scratch);
    }

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Y);
    qr_scratch.setIdentity();
    Eigen::MatrixXd Q = qr.householderQ() * qr_scratch;

    Eigen::MatrixXd AQ = apply(Q);
    Eigen::MatrixXd B  = Q.transpose() * AQ;   // k_ext x k_ext, symmetric
    AQ.resize(0, 0);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
    if (es.info() != Eigen::Success)
        throw std::runtime_error("power_iterate_and_project: eigendecomposition of projected matrix failed.");

    EighResult result;
    // es.eigenvalues() is ascending; tail(k_target).reverse() -> descending top-k.
    result.eigenvalues = es.eigenvalues().tail(k_target).reverse();
    // Materialise the reversed block into a plain contiguous matrix before the
    // n x k_target DGEMM — rowwise().reverse() on a block expression forces
    // column-by-column scatter during DGEMM, defeating tiling.
    const Eigen::MatrixXd evecs_sorted =
        es.eigenvectors().rightCols(k_target).rowwise().reverse().eval();
    result.eigenvectors = Q * evecs_sorted;
    return result;
}

// Convenience one-shot entry point: sketch + power-iterate + project.
template <typename MatVecApply>
EighResult randomized_symmetric_eigh(
    MatVecApply&& apply,
    int n,
    int k_target,
    int oversample = 20,
    int power_iter = 3,
    const Eigen::MatrixXd* warm_start = nullptr)
{
    const int k_ext = std::min(k_target + oversample, n - 1);
    auto [omega, Y] = build_randomized_sketch(apply, n, k_ext, warm_start);
    (void)omega;
    return power_iterate_and_project(std::forward<MatVecApply>(apply), std::move(Y), k_target, power_iter);
}

// ─────────────────────────────────────────────────────────────────────────
// Lanczos (Spectra), for well-separated top spectra with few matvecs
// ─────────────────────────────────────────────────────────────────────────

// Spectra requires single-vector perform_op(x_in, y_out); this wraps a
// block-capable `apply` functor (the same one rSVD uses) so call sites
// don't need to maintain two different matvec adapters for the same matrix.
template <typename MatVecApply>
struct SpectraMatVecOp {
    using Scalar = double;
    MatVecApply apply;
    int n;
    int rows() const { return n; }
    int cols() const { return n; }
    void perform_op(const double* x_in, double* y_out) const {
        Eigen::Map<const Eigen::VectorXd> x(x_in, n);
        Eigen::Map<Eigen::VectorXd>       y(y_out, n);
        y.noalias() = apply(x);
    }
};

template <typename MatVecApply>
EighResult lanczos_symmetric_eigh(
    MatVecApply&& apply,
    int n,
    int k_target,
    int ncv = -1)
{
    if (ncv <= 0) ncv = std::min(n, std::max(3 * k_target + 1, 30));

    SpectraMatVecOp<std::decay_t<MatVecApply>> op{std::forward<MatVecApply>(apply), n};
    Spectra::SymEigsSolver<SpectraMatVecOp<std::decay_t<MatVecApply>>> eigs(op, k_target, ncv);
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestAlge);
    if (eigs.info() != Spectra::CompInfo::Successful)
        throw std::runtime_error("lanczos_symmetric_eigh: Spectra eigensolver failed.");

    EighResult result;
    result.eigenvalues = eigs.eigenvalues();
    result.eigenvectors = eigs.eigenvectors();
    return result;
}

// ─────────────────────────────────────────────────────────────────────────
// Tall-skinny thin SVD via QR reduction
// ─────────────────────────────────────────────────────────────────────────

struct ThinSVDResult {
    Eigen::MatrixXd U;                // n x k, orthonormal columns
    Eigen::VectorXd singular_values;  // descending, size k
};

// Thin SVD of an already-materialised n x k matrix Z (k << n), e.g. the
// Nystrom sketch. Eigen's BDCSVD has no LAPACKE binding (unlike
// HouseholderQR and SelfAdjointEigenSolver, which do when EIGEN_USE_LAPACKE
// is set), so calling it directly on an n x k matrix leaves it doing
// bidiagonalization/divide-and-conquer work at n scale without BLAS3
// dispatch. Reducing via a blocked QR first confines that non-BLAS3-backed
// step to the k x k factor R, where its cost is negligible regardless of
// how it's implemented, while the n-scale work (QR, and the final Q * U_R)
// goes through LAPACKE_dgeqrf/dorgqr and GEMM respectively.
inline ThinSVDResult tall_skinny_thin_svd(const Eigen::MatrixXd& Z) {
    const int n = static_cast<int>(Z.rows());
    const int k = static_cast<int>(Z.cols());

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Z);
    const Eigen::MatrixXd R = qr.matrixQR().topRows(k).triangularView<Eigen::Upper>();
    const Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n, k);

    Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU> svd_r(R);

    ThinSVDResult result;
    result.singular_values = svd_r.singularValues();
    result.U = Q * svd_r.matrixU();   // n x k GEMM
    return result;
}

} // namespace gcta_eigh