/*
 * chunked_grm_matvec.hpp
 *
 * K @ X for a symmetric n x n matrix K, computed from lower-triangular tiles
 * read on demand, without ever materializing a dense n x n K in memory.
 * Companion to symmetric_eigendecomp.hpp: this produces the same `apply`
 * signature (X -> K@X) that build_randomized_sketch / power_iterate_and_project
 * / SpectraMatVecOp already consume, so it's a drop-in alternative to
 * `K_dbl * X` — nothing in symmetric_eigendecomp.hpp needs to change.
 *
 * Tiling: rows/cols [0,n) are split into fixed-size blocks B_0..B_{m-1}
 * (same partition for rows and columns). For each pair of blocks (i, j)
 * with j <= i, the tile K[B_i, B_j] is read exactly once per apply() call:
 *
 *   - Off-diagonal (j < i): the tile is used directly for Y[B_i] and,
 *     via its transpose, reflected into Y[B_j] — this is what makes a single
 *     read of the lower triangle sufficient; the "missing" upper-triangle
 *     data is never read from disk, it's supplied algebraically.
 *   - Diagonal (j == i): the tile only has its own lower triangle valid
 *     (matching how the GRM was written — see the loader's own comment on
 *     the source symmetrization), so it's locally mirrored via
 *     selfadjointView<Lower> before use.
 *
 * This is the read-side mirror of the block-tiled construction already used
 * elsewhere (off-diagonal dgemm + diagonal dsyrk) — same tile grid, same
 * lower-triangular-only I/O, applied to a matvec instead of a build.
 */
#pragma once

#include <Eigen/Dense>
#include <functional>
#include <vector>

namespace gcta_chunked {

// Reads K[rs:re, cs:ce] for a tile at or below the diagonal (cs_block <=
// rs_block, so ce <= re always holds for j==i, and the whole tile is valid
// lower-triangular data for j<i). Must return an (re-rs) x (ce-cs) matrix.
//
// Implementation note for whoever wires this to grm_binary_io.hpp: the file
// stores float32 (row-major packed lower triangle); read as float and
// .cast<double>() once here, at the tile level, rather than widening the
// whole file up front. Tiles are small and short-lived, so the widen is
// cheap regardless of tile size, and this is the only place precision needs
// to be decided — everything downstream of this callback is already double.
using TileReader = std::function<Eigen::MatrixXd(int rs, int re, int cs, int ce)>;

struct BlockPartition {
    std::vector<int> starts;  // starts[i]..starts[i+1) is block i; starts.back() == n

    explicit BlockPartition(int n, int block_size) {
        for (int s = 0; s < n; s += block_size) starts.push_back(s);
        starts.push_back(n);
    }
    int num_blocks() const { return static_cast<int>(starts.size()) - 1; }
    int block_start(int i) const { return starts[i]; }
    int block_end(int i) const { return starts[i + 1]; }
};

// K @ X without ever holding a dense n x n K. `read_tile` is called exactly
// once per (row_block, col_block) pair with col_block <= row_block — i.e.
// once per element of the lower-triangular block grid, matching how the
// matrix is stored on disk.
inline Eigen::MatrixXd chunked_symmetric_matvec(
    const TileReader& read_tile,
    int n,
    int block_size,
    const Eigen::MatrixXd& X)
{
    const BlockPartition part(n, block_size);
    const int m = part.num_blocks();

    Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(n, X.cols());

    for (int i = 0; i < m; ++i) {
        const int rs = part.block_start(i), re = part.block_end(i);
        for (int j = 0; j <= i; ++j) {
            const int cs = part.block_start(j), ce = part.block_end(j);
            Eigen::MatrixXd tile = read_tile(rs, re, cs, ce);  // (re-rs) x (ce-cs)

            if (i == j) {
                // Diagonal block: only its own lower triangle is valid data
                // (mirrors how it was written); mirror locally before use.
                Eigen::MatrixXd tile_full = tile.selfadjointView<Eigen::Lower>();
                tile.resize(0, 0);  // fully absorbed into tile_full; not needed for the GEMM below
                Y.middleRows(rs, re - rs).noalias() += tile_full * X.middleRows(rs, re - rs);
            } else {
                // Off-diagonal: tile = K[B_i, B_j]; reflect its transpose
                // into B_j's output row range to account for K[B_j, B_i]
                // without ever reading it from disk.
                Y.middleRows(rs, re - rs).noalias() += tile * X.middleRows(cs, ce - cs);
                Y.middleRows(cs, ce - cs).noalias() += tile.transpose() * X.middleRows(rs, re - rs);
            }
        }
    }
    return Y;
}

// Diagonal of K, needed for trace(K) (used by both auto-k's lambda_plus
// bookkeeping-adjacent checks and EIG99's target mass) without reading any
// off-diagonal data at all — just the m diagonal tiles' own diagonals.
inline Eigen::VectorXd chunked_diagonal(const TileReader& read_tile, int n, int block_size) {
    const BlockPartition part(n, block_size);
    const int m = part.num_blocks();
    Eigen::VectorXd d(n);
    for (int i = 0; i < m; ++i) {
        const int rs = part.block_start(i), re = part.block_end(i);
        const Eigen::MatrixXd tile = read_tile(rs, re, rs, re);  // diagonal block only
        d.segment(rs, re - rs) = tile.diagonal();
    }
    return d;
}

// trace(K^2) = sum of squares of every entry of K. Needed by the non-EIG99
// tail-variance correction (tail_d_var), which — unlike trace(K) — genuinely
// needs every off-diagonal entry, not just the diagonal. Reuses the same
// lower-triangular tile grid: a diagonal tile's contribution is the squared
// Frobenius norm of its locally-mirrored (full) form; an off-diagonal
// tile's squared Frobenius norm is counted twice (once for K[B_i,B_j], once
// for its mirror K[B_j,B_i], which has identical squared entries). This is
// one full pass over every tile — comparable cost to one apply() call —
// done once at the end of basis construction, not per escalation round.
inline double chunked_trace_K_squared(const TileReader& read_tile, int n, int block_size) {
    const BlockPartition part(n, block_size);
    const int m = part.num_blocks();
    double total = 0.0;
    for (int i = 0; i < m; ++i) {
        const int rs = part.block_start(i), re = part.block_end(i);
        for (int j = 0; j <= i; ++j) {
            const int cs = part.block_start(j), ce = part.block_end(j);
            Eigen::MatrixXd tile = read_tile(rs, re, cs, ce);
            if (i == j) {
                Eigen::MatrixXd tile_full = tile.selfadjointView<Eigen::Lower>();
                tile.resize(0, 0);
                total += tile_full.squaredNorm();
            } else {
                total += 2.0 * tile.squaredNorm();
            }
        }
    }
    return total;
}

} // namespace gcta_chunked
