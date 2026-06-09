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
#include "MLMA_stream_common.hpp"
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
#include <array>
#include <charconv>
#include <limits>
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

template <typename T>
string to_text(T x)
{
    std::array<char, 128> buf{};
    auto [ptr, ec] = std::to_chars(buf.data(), buf.data() + buf.size(),
                                   x, std::chars_format::general,
                                   std::numeric_limits<T>::max_digits10);
    if (ec == std::errc())
        return string(buf.data(), static_cast<size_t>(ptr - buf.data()));
    return to_string(x);
}

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
        st.lambda_tail_f = static_cast<float>(lambda_tail);

        // Store as k×n (transposed) for cache-efficient GEMM in the hot MLMA loop.
        Eigen::MatrixXf Uk(k, hdr.n);
        must_read(Uk.data(),
                  static_cast<std::streamsize>((size_t)hdr.n * k * sizeof(float)));

        Eigen::VectorXf dk(k);
        must_read(dk.data(), static_cast<std::streamsize>(k * sizeof(float)));
        st.dk_f = dk;

        if (!no_adj_covar) {
            st.b.resize(hdr.x_c);
            must_read(st.b.data(),
                      static_cast<std::streamsize>(hdr.x_c * sizeof(float)));
        }

        // Read varcmp to reconstruct ck
        Eigen::VectorXf vc(hdr.num_varcmp);
        must_read(vc.data(),
                  static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));
        st.varcmp = vc;

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

    st.varcmp.resize(hdr.num_varcmp);
    must_read(st.varcmp.data(),
              static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));
    return st;
}

void writeRemlState(const string& filename, const RemlState& st, bool no_adj_covar)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open())
        LOGGER.e(0, "cannot open [" + filename + "] for writing.");

    if (st.is_woodbury) {
        struct Header {
            char    magic[4] = {'T', 'U', 'N', 'A'};
            int32_t n = 0, x_c = 0, num_varcmp = 0, num_r_indx = 0;
        } hdr;
        hdr.n = st.n;
        hdr.x_c = st.x_c;
        hdr.num_varcmp = static_cast<int32_t>(st.varcmp.size());
        hdr.num_r_indx = static_cast<int32_t>(st.varcmp.size());
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

        const int32_t k = static_cast<int32_t>(st.wb.ck_f.size());
        out.write(reinterpret_cast<const char*>(&k), sizeof(int32_t));

        const double lambda_tail = static_cast<double>(st.lambda_tail_f);
        out.write(reinterpret_cast<const char*>(&lambda_tail), sizeof(double));

        if (st.wb.Uk_f.rows() != k || st.wb.Uk_f.cols() != st.n)
            LOGGER.e(0, "invalid Woodbury REML state dimensions before save.");
        if (st.dk_f.size() != k)
            LOGGER.e(0, "invalid Woodbury eigenvalue vector before save.");

        out.write(reinterpret_cast<const char*>(st.wb.Uk_f.data()),
                  static_cast<std::streamsize>(static_cast<size_t>(st.n) * k * sizeof(float)));
        out.write(reinterpret_cast<const char*>(st.dk_f.data()),
                  static_cast<std::streamsize>(k * sizeof(float)));

        if (!no_adj_covar) {
            if (st.b.size() != st.x_c)
                LOGGER.e(0, "invalid fixed-effect vector length before save.");
            out.write(reinterpret_cast<const char*>(st.b.data()),
                      static_cast<std::streamsize>(st.x_c * sizeof(float)));
        }

        out.write(reinterpret_cast<const char*>(st.varcmp.data()),
                  static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));
    } else {
        struct Header {
            char    magic[4] = {'G', 'O', 'B', 'Y'};
            int32_t n = 0, x_c = 0, num_varcmp = 0, num_r_indx = 0;
        } hdr;
        hdr.n = st.n;
        hdr.x_c = st.x_c;
        hdr.num_varcmp = static_cast<int32_t>(st.varcmp.size());
        hdr.num_r_indx = static_cast<int32_t>(st.varcmp.size());
        out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

        Eigen::MatrixXf Vi_f;
        if (st.is_llt) {
            if (st.Vi_L_f.rows() != st.n || st.Vi_L_f.cols() != st.n)
                LOGGER.e(0, "invalid LLT REML state dimensions before save.");
            Vi_f = Eigen::MatrixXf::Identity(st.n, st.n);
            st.Vi_L_f.triangularView<Eigen::Lower>().solveInPlace(Vi_f);
            st.Vi_L_f.triangularView<Eigen::Lower>().transpose().solveInPlace(Vi_f);
        } else {
            Vi_f = st.Vi;
        }

        const size_t tri = static_cast<size_t>(st.n) * (st.n + 1) / 2;
        vector<float> packed(tri);
        size_t idx = 0;
        for (int32_t j = 0; j < st.n; ++j)
            for (int32_t i = j; i < st.n; ++i)
                packed[idx++] = Vi_f(i, j);
        out.write(reinterpret_cast<const char*>(packed.data()),
                  static_cast<std::streamsize>(tri * sizeof(float)));

        if (!no_adj_covar) {
            if (st.b.size() != st.x_c)
                LOGGER.e(0, "invalid fixed-effect vector length before save.");
            out.write(reinterpret_cast<const char*>(st.b.data()),
                      static_cast<std::streamsize>(st.x_c * sizeof(float)));
        }

        out.write(reinterpret_cast<const char*>(st.varcmp.data()),
                  static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));
    }

    out.flush();
    if (!out)
        LOGGER.e(0, "write error on [" + filename + "] — disk full or I/O failure.");
}

void write_hsq_from_ctx(const string& out_prefix, const RemlCtx& ctx)
{
    const string hsq_file = out_prefix + ".hsq";
    std::ofstream o_hsq(hsq_file);
    if (!o_hsq) LOGGER.e(0, "cannot open [" + hsq_file + "] for writing.");

    o_hsq << "Source\tVariance\tSE\n";
    const int m = static_cast<int>(ctx.varcmp.size());
    for (int i = 0; i < m; ++i) {
        const string name = (i < static_cast<int>(ctx.var_name.size()))
            ? ctx.var_name[i]
            : (i + 1 == m ? "V(e)" : ("V(G" + to_string(i + 1) + ")"));
        const double se = (i < static_cast<int>(ctx.varcmp_se.size())) ? ctx.varcmp_se[i] : 0.0;
        o_hsq << name << "\t" << to_text(ctx.varcmp[i]) << "\t" << to_text(se) << "\n";
    }

    o_hsq << "Vp\t" << to_text(ctx.Vp) << "\t" << to_text(ctx.Vp_se) << "\n";

    const int ngen = static_cast<int>(ctx.hsq.size());
    for (int i = 0; i < ngen; ++i) {
        const string hname = (i < static_cast<int>(ctx.hsq_name.size()))
            ? ctx.hsq_name[i]
            : ("V(G" + to_string(i + 1) + ")/Vp");
        const double hse = (i < static_cast<int>(ctx.hsq_se.size())) ? ctx.hsq_se[i] : 0.0;
        o_hsq << hname << "\t" << to_text(ctx.hsq[i]) << "\t" << to_text(hse) << "\n";
    }

    if (ctx.has_logL)
        o_hsq << "logL\t" << to_text(ctx.logL) << "\n";

    o_hsq << "n\t" << ctx.n << "\n";
    o_hsq.close();
    LOGGER.i(0, "Summary REML results saved to [" + hsq_file + "].");
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
    const bool has_save_reml   = options_in.find("--save-reml")   != options_in.end()
                                  && !options_in["--save-reml"].empty();

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

    if (has_load_reml && has_save_reml)
        LOGGER.e(0, "--mlma-stream does not allow --load-reml with --save-reml.");

    if (has_save_reml) {
        options["save_reml"] = options_in["--save-reml"][0];
        options_in.erase("--save-reml");
    }

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

        // --reml-woodbury [k]  (k optional; -1 = auto-k, -2 = EIG99)
        if (options_in.find("--reml-woodbury") != options_in.end()) {
            const auto& vals = options_in["--reml-woodbury"];
            if (!vals.empty() && !vals[0].empty()) {
                if (vals[0] == "EIG99")
                    options_d["woodbury_rank"] = -2.0;  // Eigenvalue mass k
                else if (vals[0] == "auto")
                    options_d["woodbury_rank"] = -1.0;  // auto-k
                else
                    options_d["woodbury_rank"] = std::stod(vals[0]);
            } else
                options_d["woodbury_rank"] = -1.0;   // auto-k
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

            // Compact phenotype values to the kept sample order.
            phenos_vec = compact_sample_vector(phenos_vec, remain_index);
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
            ctx.grm_N.resize(1, 1);
            ctx.grm_N(0, 0) = m_all;
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

            // Reuse the already-computed REML summary to write .hsq (same feature as --mlma).
            write_hsq_from_ctx(out_prefix, ctx);

            state = reml::build_reml_state(ctx);

            if (options.count("save_reml")) {
                const string save_reml_file = options.at("save_reml");
                LOGGER.i(0, "Saving REML state to [" + save_reml_file + "] ...");
                writeRemlState(save_reml_file, state, no_adj_covar);
                LOGGER.i(0, "REML estimation completed. Use --load-reml to perform association tests.");
                delete geno;
                delete marker;
                delete pheno;
                continue;
            }

            // y_adj = y - X*b
            {
                const Eigen::MatrixXf Xf = X_design.cast<float>();
                y_vec -= Xf * state.b;
            }
        }

        // ---- Open output file ----
        const string out_file = out_prefix + ".mlma";
        std::ofstream ofile(out_file);
        if (!ofile) LOGGER.e(0, "cannot open [" + out_file + "] for writing.");
        ofile << "Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\t"
              << (log_pval ? "log_p" : "p") << "\n";

        // ---- Stream SNPs ----
        run_mlma_stream_association(state, y_vec, geno, marker, n, log_pval, ofile);
        LOGGER << "\nAssociation results saved to [" << out_file << "]." << std::endl;
        ofile.close();

        delete geno;
        delete marker;
        delete pheno;
    }
}
