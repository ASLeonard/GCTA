#pragma once

#include "Geno.h"
#include "Marker.h"
#include "Logger.h"
#include "RemlState.hpp"
#include "cpu.h"
#include "mlma_woodbury.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <type_traits>
#include <span>
#include <string>
#include <vector>

template <typename T, typename IndexT,
          typename = std::enable_if_t<std::is_integral_v<IndexT>>>
inline std::vector<T> compact_sample_vector(const std::vector<T>& values,
                                            const std::vector<IndexT>& keep)
{
    std::vector<T> compacted;
    compacted.reserve(keep.size());
    for (IndexT index : keep) compacted.push_back(values[static_cast<size_t>(index)]);
    return compacted;
}

template <typename IndexT,
          typename = std::enable_if_t<std::is_integral_v<IndexT>>>
inline Eigen::VectorXd compact_sample_vector(const Eigen::VectorXd& values,
                                             const std::vector<IndexT>& keep)
{
    Eigen::VectorXd compacted(keep.size());
    for (size_t i = 0; i < keep.size(); ++i)
        compacted[static_cast<Eigen::Index>(i)] = values[static_cast<Eigen::Index>(keep[i])];
    return compacted;
}

template <typename IndexT,
          typename = std::enable_if_t<std::is_integral_v<IndexT>>>
inline Eigen::MatrixXd compact_sample_rows(const Eigen::MatrixXd& values,
                                           const std::vector<IndexT>& keep)
{
    Eigen::MatrixXd compacted(keep.size(), values.cols());
    for (size_t i = 0; i < keep.size(); ++i)
        compacted.row(static_cast<Eigen::Index>(i)) = values.row(static_cast<Eigen::Index>(keep[i]));
    return compacted;
}

inline void run_mlma_stream_association(RemlState& state,
                                        const Eigen::VectorXf& y_adj,
                                        const Eigen::VectorXf& w_sqrt,
                                        Geno* geno,
                                        Marker* marker,
                                        int n,
                                        bool log_pval,
                                        std::ofstream& ofile)
{
    Eigen::VectorXf Vi_y(n);
    Eigen::LLT<Eigen::MatrixXf> Vi_llt;

    const bool use_wb  = state.is_woodbury;
    const bool use_llt = state.is_llt;

    if (use_wb) {
        Vi_y = woodbury_apply_Vi_f(state.wb, y_adj);
    } else if (use_llt) {
        // V^{-1} y = L^{-T}(L^{-1} y) via two float triangular solves.
        // Avoids dpotri + a second Cholesky decomposition (saves 2×O(n³/3)).
        Vi_y = y_adj;
        cblas_strsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                    n, state.Vi_L_f.data(), n, Vi_y.data(), 1);
        cblas_strsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
                    n, state.Vi_L_f.data(), n, Vi_y.data(), 1);
    } else {
        Vi_y.noalias() = state.Vi * y_adj;
        Vi_llt = Eigen::LLT<Eigen::MatrixXf>(state.Vi);
        if (Vi_llt.info() != Eigen::Success)
            LOGGER.e(0, "REML state Vi matrix is not positive definite.");
        state.Vi.resize(0, 0);  // free RAM immediately
    }

    // Apply sqrt(W) to Vi_y once; genotype columns are scaled per-SNP in the callback.
    Vi_y.array() *= w_sqrt.array();

    const uint32_t total_m = marker->count_extract();
    LOGGER << "\nRunning association tests for " << total_m << " SNPs..." << std::endl;

    constexpr int BLOCK = 10000;
    Eigen::MatrixXf X_block(n, BLOCK);
    Eigen::VectorXf Xt_Vi_y(BLOCK);
    Eigen::VectorXf xvx_diag(BLOCK);

    int  snp_done = 0;
    int  last_pct = -1;
    const WoodburyMLMACache& wb = state.wb;

    // Pre-allocate per-block buffers (reused across callbacks, avoids heap churn).
    // GenoBufItem: item.geno is resized to n on the first getGenoDouble call and
    // reused thereafter — same pattern as GRM.cpp's gbufitems array.
    std::vector<GenoBufItem> gbuf_items(BLOCK);
    std::vector<uint8_t>     valid_v(BLOCK, 0);  // uint8_t avoids vector<bool> bit-packing
    std::vector<float>       af_v(BLOCK, 0.0f);

    // Batch output buffer to reduce write() syscall overhead.
    std::vector<char> io_buf(4 << 20);  // 4 MiB
    size_t io_pos = 0;
    std::vector<char> line_scratch;
    auto flush_io = [&]() {
        if (io_pos > 0) {
            ofile.write(io_buf.data(), static_cast<std::streamsize>(io_pos));
            io_pos = 0;
        }
    };
    auto append_io = [&](const char* data, size_t len) {
        if (len == 0) return;
        if (len > io_buf.size()) {
            flush_io();
            ofile.write(data, static_cast<std::streamsize>(len));
            return;
        }
        if (io_pos + len > io_buf.size()) flush_io();
        std::memcpy(io_buf.data() + io_pos, data, len);
        io_pos += len;
    };

    auto callback = [&](uintptr_t* buf, std::span<const uint32_t> exIdx) {
        const int bs = static_cast<int>(exIdx.size());

        for (int i = 0; i < bs; ++i) {
            valid_v[i] = 0;
            af_v[i]    = 0.0f;
            GenoBufItem& item = gbuf_items[i];
            item.extractedMarkerIndex = exIdx[i];  // MUST be set before getGenoDouble
            geno->getGenoDouble(buf, i, &item);
            if (!item.valid) { X_block.col(i).setZero(); continue; }
            if (static_cast<int>(item.geno.size()) != n) {
                LOGGER.e(0, "internal error: SNP " + std::to_string(exIdx[i])
                            + " returned geno.size=" + std::to_string(item.geno.size())
                            + " but expected " + std::to_string(n) + ".");
            }
            valid_v[i] = 1;
            af_v[i]    = static_cast<float>(item.af);
            X_block.col(i) =
                Eigen::Map<const Eigen::VectorXd>(item.geno.data(), n).cast<float>();
            X_block.col(i).array() *= w_sqrt.array();
        }


        // Xt_Vi_y = X^T Vi y  (computed BEFORE STRMM overwrites X_block)
        Xt_Vi_y.head(bs).noalias() =
            X_block.leftCols(bs).transpose() * Vi_y;

        // xvx_diag = diag(X^T Vi^{-1} X)
        if (use_wb) {
            woodbury_xvx_diag_block(wb, X_block, bs, xvx_diag);
        } else if (use_llt) {
            // x^T V^{-1} x = ||L^{-1} x||^2 (L = lower Cholesky of V).
            // In-place STRSM overwrites X_block with L^{-1} X_block.
            cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                        n, bs, 1.0f,
                        state.Vi_L_f.data(), n,
                        X_block.data(), n);
            xvx_diag.head(bs) = X_block.leftCols(bs).colwise().squaredNorm();
        } else {
            // x^T V^{-1} x = ||L_vi^T x||^2 (L_vi = Cholesky of V^{-1}).
            cblas_strmm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit,
                        n, bs, 1.0f,
                        Vi_llt.matrixLLT().data(), n,
                        X_block.data(), n);
            xvx_diag.head(bs) = X_block.leftCols(bs).colwise().squaredNorm();
        }


        // Write per-SNP results into batch output buffer to minimise write() overhead.
        for (int i = 0; i < bs; ++i) {
            const uint32_t raw      = marker->getRawIndex(exIdx[i]);
            const std::string& chr  = marker->getRawChr(raw);
            const std::string& name = marker->getRawName(raw);
            const unsigned     bp   = static_cast<unsigned>(marker->getRawBp(raw));
            const bool eff_rev      = marker->isEffecRev(exIdx[i]);
            const std::string& a1   = eff_rev ? marker->getRawA2(raw) : marker->getRawA1(raw);
            const std::string& a2   = eff_rev ? marker->getRawA1(raw) : marker->getRawA2(raw);

            float beta_val = 0.0f, se_val = 0.0f;
            double pval_val = 0.0;
            const bool stat_ok = valid_v[i] &&
                mlma_snp_stat(Xt_Vi_y[i], xvx_diag[i], log_pval,
                            beta_val, se_val, pval_val);

            if (!stat_ok) {
                std::format_to(std::back_inserter(io_buf),
                    "{}\t{}\t{}\t{}\t{}\tNA\tNA\tNA\tNA\n",
                    chr, name, bp, a1, a2);
            } else {
                std::format_to(std::back_inserter(io_buf),
                    "{}\t{}\t{}\t{}\t{}\t{:.6g}\t{:.6g}\t{:.6g}\t{:.6g}\n",
                    chr, name, bp, a1, a2,
                    static_cast<double>(af_v[i]),
                    static_cast<double>(beta_val),
                    static_cast<double>(se_val),
                    pval_val);
            }
        }

        snp_done += bs;
        const int cur_pct = (total_m > 0)
            ? static_cast<int>((uint64_t)snp_done * 100 / total_m) : 100;
        if (cur_pct != last_pct) {
            LOGGER.p(0, std::to_string(snp_done) + " / " + std::to_string(total_m)
                        + " SNPs (" + std::to_string(cur_pct) + "%)");
            last_pct = cur_pct;
        }
    };

    const std::vector<uint32_t>& extractIndex = marker->get_extract_index();
    geno->loopDouble(extractIndex, BLOCK,
                     /*bMakeGeno*/   true,
                     /*bGenoCenter*/ true,
                     /*bGenoStd*/    false,
                     /*bMakeMiss*/   true,
                     {callback});

    flush_io();
}