/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Streaming MLMA (--MLMA) – loads REML state from a binary file produced by
 * "--mlma --save-reml" and streams SNPs via Geno::loopDouble instead of
 * pre-loading the full genotype matrix into RAM.
 *
 * Mandatory flag: --load-reml <file>
 */

#include "MLMA.h"
#include "Pheno.h"
#include "Marker.h"
#include "Geno.h"
#include "Covar.h"
#include "Logger.h"
#include "cpu.h"          // cblas_strmm / cblas_strsv / cblas_strsm
#include "main/StatFunc.h"
#include "mlma_woodbury.hpp"
#include "RemlState.hpp"
#include "RemlCtx.hpp"
#include "RemlEngine.hpp"
#include "grm_binary_io.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <cstring>

using std::map;
using std::string;
using std::to_string;
using std::vector;
using std::function;

// ---------------------------------------------------------------------------
// Static member definitions
// ---------------------------------------------------------------------------
map<string, string>  MLMA::options;
map<string, double>  MLMA::options_d;
map<string, vector<double>> MLMA::options_vd;
vector<string>       MLMA::processFunctions;

// ---------------------------------------------------------------------------
// Internal types (anonymous namespace)
// ---------------------------------------------------------------------------
namespace {

// Shared GRM I/O helpers from include/grm_binary_io.hpp.
using gcta_grm_io::read_grm_binary;
using gcta_grm_io::match_ids_to_grm;

// RemlState is now defined in include/RemlState.hpp (shared with MLMA_loco).
// The readRemlState() function below remains file-local (only MLMA_stream needs it).

// Read the binary REML state written by save_reml_state().
// When !no_adj_covar the 'b' vector is loaded; otherwise it is skipped.
RemlState readRemlState(const string& filename, bool no_adj_covar)
{
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
        LOGGER.e(0, "cannot open [" + filename + "] — run --mlma --save-reml first.");

    auto must_read = [&](void* dst, std::streamsize nbytes) {
        f.read(reinterpret_cast<char*>(dst), nbytes);
        if (f.gcount() != nbytes)
            LOGGER.e(0, "unexpected EOF in REML state file [" + filename + "].");
    };

    struct Header {
        char    magic[4];
        int32_t n, x_c, num_varcmp, num_r_indx;
    } hdr;
    must_read(&hdr, sizeof(hdr));

    if (hdr.n <= 0 || hdr.x_c < 0 || hdr.num_varcmp <= 0)
        LOGGER.e(0, "[" + filename + "] has invalid header dimensions.");

    const std::string_view magic(hdr.magic, 4);
    RemlState st;
    st.n   = hdr.n;
    st.x_c = hdr.x_c;

    // ------------------------------------------------------------------ TUNA
    if (magic == "TUNA") {
        st.is_woodbury = true;
        int32_t k = 0;
        must_read(&k, sizeof(int32_t));
        if (k <= 0)
            LOGGER.e(0, "[" + filename + "] has invalid Woodbury rank k=" + to_string(k) + ".");

        double lambda_tail = 0.0;
        must_read(&lambda_tail, sizeof(double));

        // Store as k×n (transposed) for cache-efficient GEMM in the hot MLMA loop.
        Eigen::MatrixXf Uk(k, hdr.n);
        must_read(Uk.data(),
                  static_cast<std::streamsize>((size_t)hdr.n * k * sizeof(float)));

        Eigen::VectorXf dk(k);
        must_read(dk.data(), static_cast<std::streamsize>(k * sizeof(float)));

        if (!no_adj_covar) {
            st.b.resize(hdr.x_c);
            must_read(st.b.data(),
                      static_cast<std::streamsize>(hdr.x_c * sizeof(float)));
        }

        // Read varcmp to reconstruct ck
        Eigen::VectorXf vc(hdr.num_varcmp);
        must_read(vc.data(),
                  static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));

        const double sg2    = static_cast<double>(vc[0]);
        const double se2    = static_cast<double>(vc[hdr.num_varcmp - 1]);
        const double sig2e  = sg2 * lambda_tail + se2;

        st.wb.Uk_f.swap(Uk);
        st.wb.ck_f.resize(k);
        for (int j = 0; j < k; ++j) {
            double delta   = std::max(0.0, static_cast<double>(dk[j]) - lambda_tail);
            double sd      = sg2 * delta;
            st.wb.ck_f[j]  = static_cast<float>(sd / (sig2e + sd));
        }
        st.wb.sqrt_ck_f    = st.wb.ck_f.cwiseSqrt();
        st.wb.sigma2_eff_f = static_cast<float>(sig2e);
        return st;
    }

    // ------------------------------------------------------------------ GOBY
    if (magic != "GOBY")
        LOGGER.e(0, "[" + filename + "] is not a valid REML state file (bad magic bytes).");

    const size_t tri = static_cast<size_t>(hdr.n) * (hdr.n + 1) / 2;
    vector<float> buf(tri);
    must_read(buf.data(), static_cast<std::streamsize>(tri * sizeof(float)));

    // Unpack lower-triangle into full symmetric matrix
    st.Vi.resize(hdr.n, hdr.n);
    size_t idx = 0;
    for (int32_t j = 0; j < hdr.n; ++j)
        for (int32_t i = j; i < hdr.n; ++i, ++idx) {
            st.Vi(i, j) = buf[idx];
            st.Vi(j, i) = buf[idx];
        }

    if (!no_adj_covar) {
        st.b.resize(hdr.x_c);
        must_read(st.b.data(),
                  static_cast<std::streamsize>(hdr.x_c * sizeof(float)));
    }
    return st;
}

} // anonymous namespace


// ---------------------------------------------------------------------------
// registerOption
// ---------------------------------------------------------------------------
int MLMA::registerOption(map<string, vector<string>>& options_in)
{
    // Always capture --out (shared flag)
    // Use the pre-processed "out" key (no dashes), which is always set by main.cpp
    if (!options_in["out"].empty())
        options["out"] = options_in["out"][0];

    const bool has_mlma_stream = options_in.find("--mlma-stream") != options_in.end();
    const bool has_load_reml   = options_in.find("--load-reml")   != options_in.end()
                                  && !options_in["--load-reml"].empty();

    auto capture_common_reml_flags = [&]() {
        if (options_in.find("--reml-alg") != options_in.end()
                && !options_in["--reml-alg"].empty()) {
            options_d["reml_alg"] = std::stod(options_in["--reml-alg"][0]);
            options_in.erase("--reml-alg");
        }
        if (options_in.find("--reml-no-constrain") != options_in.end()) {
            options["no_constrain"] = "1";
            options_in.erase("--reml-no-constrain");
        }
        if (options_in.find("--reml-diagV-adj") != options_in.end()
                && !options_in["--reml-diagV-adj"].empty()) {
            options_d["reml_diagV_adj"] = std::stod(options_in["--reml-diagV-adj"][0]);
            options_in.erase("--reml-diagV-adj");
        }
    };

    if (!has_mlma_stream)
        return 0;

    capture_common_reml_flags();

    if (has_load_reml) {
        // Pre-built REML state: load from file.
        options["load_reml"] = options_in["--load-reml"][0];
        options_in.erase("--load-reml");
        // REML tuning flags are irrelevant when a saved state is loaded; consume to suppress warnings.
        options_in.erase("--reml-woodbury");
        options_in.erase("--reml-trace-approx");
        options_in.erase("--reml-maxit");
        options_in.erase("--reml-priors");
        options_in.erase("--reml-priors-var");
    } else {
        // Inline REML path: --grm is required.
        const bool has_grm = options_in.find("--grm") != options_in.end()
                              && !options_in["--grm"].empty();
        if (!has_grm)
            LOGGER.e(0, "--mlma-stream requires either --load-reml <file> or --grm <prefix>.");
        options["grm"] = options_in["--grm"][0];
        options_in.erase("--grm");

        // Optional REML tuning flags
        if (options_in.find("--reml-woodbury") != options_in.end()) {
            const auto& v = options_in["--reml-woodbury"];
            options_d["woodbury_rank"] = (!v.empty() && !v[0].empty()) ? std::stod(v[0]) : -1.0;
            options_in.erase("--reml-woodbury");
        }
        if (options_in.find("--reml-trace-approx") != options_in.end()) {
            options["trace_approx"] = "1";
            const auto& v = options_in["--reml-trace-approx"];
            if (!v.empty() && !v[0].empty())
                options_d["trace_approx_nprobes"] = std::stod(v[0]);
            options_in.erase("--reml-trace-approx");
        }
        if (options_in.find("--reml-maxit") != options_in.end()
                && !options_in["--reml-maxit"].empty()) {
            options_d["reml_maxit"] = std::stod(options_in["--reml-maxit"][0]);
            options_in.erase("--reml-maxit");
        }
        if (options_in.find("--reml-priors-var") != options_in.end()) {
            options["reml_priors_var"] = "1";
            const auto& vals = options_in["--reml-priors-var"];
            for (const auto& s : vals) {
                if (!s.empty())
                    options_vd["reml_priors_var"].push_back(std::stod(s));
            }
            options_in.erase("--reml-priors-var");
        }
        if (options_in.find("--reml-priors") != options_in.end()) {
            const auto& vals = options_in["--reml-priors"];
            for (const auto& s : vals) {
                if (!s.empty())
                    options_vd["reml_priors"].push_back(std::stod(s));
            }
            options_in.erase("--reml-priors");
        }
    }

    if (options_in.find("--mlma-no-preadj-covar") != options_in.end()) {
        options["no_adj_covar"] = "1";
        options_in.erase("--mlma-no-preadj-covar");
    }
    if (options_in.find("--log-pval") != options_in.end()) {
        options["log_pval"] = "1";
        options_in.erase("--log-pval");
    }

    processFunctions.push_back("MLMA");
    options_in.erase("--mlma-stream");
    return 1;
}


// ---------------------------------------------------------------------------
// processMain
// ---------------------------------------------------------------------------
void MLMA::processMain()
{
    for (const auto& pf : processFunctions) {
        if (pf != "MLMA") continue;

        const string out_prefix     = options.at("out");
        const bool   no_adj_covar   = options.count("no_adj_covar") > 0;
        const bool   log_pval       = options.count("log_pval") > 0;

        if (no_adj_covar)
            LOGGER.e(0, "--mlma-no-preadj-covar is not yet supported in --MLMA. "
                        "Re-run without this flag (pre-adjustment is the default).");

        // ---- Pheno / Marker / Geno (Geno re-reads pheno state in loopDouble) ----
        Pheno*  pheno  = new Pheno();
        Marker* marker = new Marker();
        Geno*   geno   = new Geno(pheno, marker);

        // ---- Phenotype ----
        const uint32_t n_init = pheno->count_keep();
        vector<string> ids;
        vector<double> phenos_vec;
        pheno->get_pheno(ids, phenos_vec);
        if (ids.size() != n_init)
            LOGGER.e(0, "--MLMA requires --pheno to be specified.");

        // ---- Covariates ----
        Covar covar;
        bool  has_covar = false;
        int   x_c       = 1;
        Eigen::MatrixXd X_design;

        {
            vector<uint32_t> remain_index, covar_index;
            has_covar = covar.getCommonSampleIndex(ids, remain_index, covar_index);

            if (has_covar) {
                LOGGER.i(0, to_string(remain_index.size()) +
                            " overlapping individuals with non-missing covariate data.");
                pheno->filter_keep_index(remain_index);
            } else if (covar.hasCovar()) {
                LOGGER.e(0, "no overlapping individuals with non-missing covariate data.");
            } else {
                // No covariate files specified – keep all phenotype samples
                remain_index.resize(ids.size());
                std::iota(remain_index.begin(), remain_index.end(), 0u);
            }

            // Reorder phenotype vector to match filtered index_keep order
            const int n_new = static_cast<int>(remain_index.size());
            vector<double> phenos2(n_new);
            for (int i = 0; i < n_new; ++i)
                phenos2[i] = phenos_vec[remain_index[i]];
            phenos_vec = std::move(phenos2);
        }

        const int n = static_cast<int>(pheno->count_keep());
        LOGGER.i(0, "After matching all files, " + to_string(n) +
                    " individuals to be included in --MLMA.");

        // Build design matrix X  (n × x_c, column-major)
        if (has_covar) {
            vector<string> final_ids = pheno->get_id(0, n - 1, "\t");
            vector<double> covar_X;
            vector<uint32_t> ck_dummy;
            covar.getCovarX(final_ids, covar_X, ck_dummy);
            const int covar_cols = (n > 0) ? static_cast<int>(covar_X.size()) / n : 0;
            x_c = 1 + covar_cols;
            X_design.resize(n, x_c);
            X_design.col(0).setOnes();
            if (covar_cols > 0)
                X_design.rightCols(covar_cols) =
                    Eigen::Map<const Eigen::MatrixXd>(covar_X.data(), n, covar_cols);
        } else {
            X_design.resize(n, 1);
            X_design.setOnes();
        }

        // ---- Obtain REML state (either from file or inline) ----
        RemlState state;
        Eigen::VectorXf y_vec(n);
        for (int i = 0; i < n; ++i) y_vec[i] = static_cast<float>(phenos_vec[i]);

        if (options.count("load_reml")) {
            const string load_reml_file = options.at("load_reml");
            LOGGER.i(0, "Loading REML state from [" + load_reml_file + "]...");
            state = readRemlState(load_reml_file, no_adj_covar);

            if (state.n != n)
                LOGGER.e(0, "Sample size mismatch: REML state n=" + to_string(state.n) +
                            " vs dataset n=" + to_string(n) +
                            ". Use the same filters (--keep/--remove/--pheno) in both runs.");
            if (state.x_c != x_c)
                LOGGER.e(0, "Covariate count mismatch: REML state x_c=" + to_string(state.x_c) +
                            " vs current x_c=" + to_string(x_c) +
                            ". Use the same --qcovar/--covar as during --save-reml.");

            // y_adj = y - X*b
            {
                const Eigen::MatrixXf Xf = X_design.cast<float>();
                y_vec -= Xf * state.b;
            }
        } else {
            // ---- Inline REML: read GRM, run reml::compute(), build RemlState ----
            const string grm_pfx = options.at("grm");
            LOGGER.i(0, "Running inline REML using GRM [" + grm_pfx + "] ...");

            vector<string> grm_ids;
            Eigen::MatrixXd G_n;
            double m_all = 0.0;
            read_grm_binary(grm_pfx, grm_ids, G_n, m_all);

            // Get post-filter analysis IDs (FID\tIID) and match to GRM
            const vector<string> analysis_ids = pheno->get_id(0, n - 1, "\t");
            const vector<int>    kp           = match_ids_to_grm(analysis_ids, grm_ids);
            for (int i = 0; i < n; ++i)
                if (kp[i] < 0)
                    LOGGER.e(0, "Individual [" + analysis_ids[i] +
                                "] not found in GRM [" + grm_pfx + "]. "
                                "Re-build the GRM from the same sample set.");

            // Subset GRM to the n analysis individuals.
            // Fast path: if GRM sample == analysis sample (in order), no copy needed.
            {
                const int  n_grm        = static_cast<int>(grm_ids.size());
                bool       is_identity  = (n_grm == n);
                for (int i = 0; i < n && is_identity; ++i)
                    if (kp[i] != i) is_identity = false;

                if (!is_identity) {
                    Eigen::MatrixXd G_sub(n, n);
                    // Column-first traversal for column-major Eigen storage.
                    for (int j = 0; j < n; ++j) {
                        const int src_col = kp[j];
                        for (int i = 0; i < n; ++i)
                            G_sub(i, j) = G_n(kp[i], src_col);
                    }
                    G_n = std::move(G_sub);
                }
                // else: G_n already contains the right n×n block
            }

            // REML tuning parameters
            const int  woodbury_rank  = options_d.count("woodbury_rank")
                ? static_cast<int>(options_d.at("woodbury_rank")) : 0;
            const bool trace_approx   = options.count("trace_approx") > 0;
            const int  trace_nprobes  = options_d.count("trace_approx_nprobes")
                ? static_cast<int>(options_d.at("trace_approx_nprobes")) : 100;
            const int  reml_maxit     = options_d.count("reml_maxit")
                ? static_cast<int>(options_d.at("reml_maxit")) : 100;
            const int  reml_alg       = options_d.count("reml_alg")
                ? static_cast<int>(options_d.at("reml_alg")) : 0;
            const int  reml_diagV_adj = options_d.count("reml_diagV_adj")
                ? static_cast<int>(options_d.at("reml_diagV_adj")) : 0;
            const bool no_constrain   = options.count("no_constrain") > 0;

            if (reml_alg < 0 || reml_alg > 2)
                LOGGER.e(0, "--reml-alg should be 0, 1 or 2.");
            if (reml_diagV_adj < 0 || reml_diagV_adj > 2)
                LOGGER.e(0, "--reml-diagV-adj should be 0, 1, or 2.");
            if (woodbury_rank != 0 && reml_alg == 1)
                LOGGER.e(0, "--reml-woodbury is incompatible with Fisher-scoring REML (--reml-alg 1). Use AI-REML (default) or EM-REML (--reml-alg 2).");

            const vector<double> priors =
                options_vd.count("reml_priors") ? options_vd.at("reml_priors") : vector<double>{};
            const vector<double> priors_var =
                options_vd.count("reml_priors_var") ? options_vd.at("reml_priors_var") : vector<double>{};

            // Build REML context
            RemlCtx ctx;
            ctx.n   = n;
            ctx.X_c = x_c;
            ctx.X   = X_design;
            {
                Eigen::VectorXd y_d(n);
                for (int i = 0; i < n; ++i) y_d[i] = phenos_vec[i];
                ctx.y     = y_d;
                ctx.y_Ssq = ctx.y.squaredNorm();
            }
            ctx.A.resize(2);
            ctx.A[0] = std::move(G_n);
            // ctx.A[1] left default (size 0 == identity convention for residual)
            ctx.r_indx   = {0, 1};
            ctx.var_name = {"V(G)", "V(e)"};
            ctx.hsq_name = {"V(G)/Vp"};
            ctx.out      = out_prefix;

            ctx.reml_mtd                 = reml_alg;
            ctx.reml_max_iter            = reml_maxit;
            ctx.reml_inv_mtd             = 0;  // LLT
            ctx.reml_diagV_adj           = reml_diagV_adj;
            ctx.woodbury_rank            = woodbury_rank;
            ctx.woodbury_buffer_factor   = 1.5;
            ctx.reml_trace_approx        = trace_approx;
            ctx.reml_trace_approx_nprobes = trace_nprobes;

            reml::compute(ctx, priors, priors_var, no_constrain);

            LOGGER.i(0, "Inline REML complete. Variance components:");
            for (size_t ci = 0; ci < ctx.varcmp.size(); ++ci)
                LOGGER.i(0, "  " + ctx.var_name[ci] + " = " + to_string(ctx.varcmp[ci]));

            state = reml::build_reml_state(ctx);

            // y_adj = y - X*b
            {
                const Eigen::MatrixXf Xf = X_design.cast<float>();
                y_vec -= Xf * state.b;
            }
        }

        // ---- Pre-scan setup ----
        Eigen::VectorXf Vi_y(n);
        Eigen::LLT<Eigen::MatrixXf> Vi_llt;

        const bool use_wb  = state.is_woodbury;
        const bool use_llt = state.is_llt;

        if (use_wb) {
            Vi_y = woodbury_apply_Vi_f(state.wb, y_vec);
        } else if (use_llt) {
            // V^{-1} y = L^{-T}(L^{-1} y) via two float triangular solves.
            // Avoids dpotri + a second Cholesky decomposition (saves 2×O(n³/3)).
            Vi_y = y_vec;
            cblas_strsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                        n, state.Vi_L_f.data(), n, Vi_y.data(), 1);
            cblas_strsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
                        n, state.Vi_L_f.data(), n, Vi_y.data(), 1);
        } else {
            Vi_y.noalias() = state.Vi * y_vec;
            Vi_llt = Eigen::LLT<Eigen::MatrixXf>(state.Vi);
            if (Vi_llt.info() != Eigen::Success)
                LOGGER.e(0, "REML state Vi matrix is not positive definite.");
            state.Vi.resize(0, 0);  // free RAM immediately
        }

        // ---- Open output file ----
        const string out_file = out_prefix + ".mlma";
        std::ofstream ofile(out_file);
        if (!ofile) LOGGER.e(0, "cannot open [" + out_file + "] for writing.");
        ofile << "Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\t"
              << (log_pval ? "log_p" : "p") << "\n";

        // ---- Stream SNPs ----
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
        vector<GenoBufItem> gbuf_items(BLOCK);
        vector<uint8_t>     valid_v(BLOCK, 0);  // uint8_t avoids vector<bool> bit-packing
        vector<float>       af_v(BLOCK, 0.0f);

        // Woodbury temp is allocated inside woodbury_xvx_diag_block; no pre-alloc needed.
        // Batch output buffer to reduce write() syscall overhead.
        std::vector<char> io_buf(4 << 20);  // 4 MiB
        size_t io_pos = 0;
        auto flush_io = [&]() {
            if (io_pos > 0) {
                ofile.write(io_buf.data(), static_cast<std::streamsize>(io_pos));
                io_pos = 0;
            }
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
                valid_v[i] = 1;
                af_v[i]    = static_cast<float>(item.af);
                X_block.col(i) =
                    Eigen::Map<const Eigen::VectorXd>(item.geno.data(), n).cast<float>();
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
                xvx_diag.head(bs) =
                    X_block.leftCols(bs).colwise().squaredNorm();
            } else {
                // Dense: in-place STRMM, then ||L^T x||^2
                cblas_strmm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasNonUnit,
                            n, bs, 1.0f,
                            Vi_llt.matrixLLT().data(), n,
                            X_block.data(), n);
                xvx_diag.head(bs) =
                    X_block.leftCols(bs).colwise().squaredNorm();
            }

            // Write per-SNP results into batch output buffer to minimise write() overhead.
            for (int i = 0; i < bs; ++i) {
                const uint32_t raw = marker->getRawIndex(exIdx[i]);
                const std::string& chr  = marker->getRawChr(raw);
                const std::string& name = marker->getRawName(raw);
                const unsigned     bp   = static_cast<unsigned>(marker->getRawBp(raw));
                const std::string& a1   = marker->getRawA1(raw);
                const std::string& a2   = marker->getRawA2(raw);

                char line[512];
                int len;
                float beta_val, se_val;
                double pval_val;
                if (!valid_v[i] || !mlma_snp_stat(Xt_Vi_y[i], xvx_diag[i], log_pval,
                                                   beta_val, se_val, pval_val)) {
                    len = std::snprintf(line, sizeof(line),
                                        "%s\t%s\t%u\t%s\t%s\tNA\tNA\tNA\tNA\n",
                                        chr.c_str(), name.c_str(), bp,
                                        a1.c_str(), a2.c_str());
                } else {
                    len = std::snprintf(line, sizeof(line),
                                        "%s\t%s\t%u\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\n",
                                        chr.c_str(), name.c_str(), bp,
                                        a1.c_str(), a2.c_str(),
                                        static_cast<double>(af_v[i]),
                                        static_cast<double>(beta_val),
                                        static_cast<double>(se_val),
                                        pval_val);
                }
                if (static_cast<size_t>(len) + io_pos > io_buf.size()) flush_io();
                std::memcpy(io_buf.data() + io_pos, line, static_cast<size_t>(len));
                io_pos += static_cast<size_t>(len);
            }

            snp_done += bs;
            const int cur_pct = (total_m > 0)
                ? static_cast<int>((uint64_t)snp_done * 100 / total_m) : 100;
            if (cur_pct != last_pct) {
                LOGGER.p(0, to_string(snp_done) + " / " + to_string(total_m)
                            + " SNPs (" + to_string(cur_pct) + "%)");
                last_pct = cur_pct;
            }
        };

        const vector<uint32_t>& extractIndex = marker->get_extract_index();
        geno->loopDouble(extractIndex, BLOCK,
                         /*bMakeGeno*/   true,
                         /*bGenoCenter*/ true,
                         /*bGenoStd*/    false,
                         /*bMakeMiss*/   true,
                         {callback});

        flush_io();  // flush remaining buffered output
        LOGGER << "\nAssociation results saved to [" << out_file << "]." << std::endl;
        ofile.close();

        delete geno;
        delete marker;
        delete pheno;
    }
}
