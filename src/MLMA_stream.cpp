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
#include "cpu.h"          // cblas_strmm
#include "main/StatFunc.h"
#include "mlma_woodbury.hpp"

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
vector<string>       MLMA::processFunctions;

// ---------------------------------------------------------------------------
// Internal types (anonymous namespace)
// ---------------------------------------------------------------------------
namespace {

struct RemlState {
    int32_t n            = 0;
    int32_t x_c          = 0;
    bool    is_woodbury  = false;
    Eigen::MatrixXf Vi;   // GOBY: full n×n V^{-1}
    WoodburyMLMACache wb; // TUNA: Woodbury low-rank factors
    Eigen::VectorXf b;    // fixed-effect coefficients (x_c elements)
};

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
    const bool has_load_reml   = options_in.find("--load-reml")   != options_in.end();

    if (!has_mlma_stream)
        return 0;

    // --load-reml is mandatory
    if (!has_load_reml || options_in["--load-reml"].empty())
        LOGGER.e(0, "--mlma-stream requires --load-reml <file>.");
    options["load_reml"] = options_in["--load-reml"][0];
    options_in.erase("--load-reml");

    if (options_in.find("--mlma-no-preadj-covar") != options_in.end()) {
        options["no_adj_covar"] = "1";
        options_in.erase("--mlma-no-preadj-covar");
    }
    if (options_in.find("--log-pval") != options_in.end()) {
        options["log_pval"] = "1";
        options_in.erase("--log-pval");
    }

    // REML flags irrelevant once --load-reml is used; consume so they don't confuse V2.
    options_in.erase("--reml-trace-approx");

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

        const string load_reml_file = options.at("load_reml");
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

        // ---- Load REML state ----
        LOGGER.i(0, "Loading REML state from [" + load_reml_file + "]...");
        RemlState state = readRemlState(load_reml_file, no_adj_covar);

        if (state.n != n)
            LOGGER.e(0, "Sample size mismatch: REML state n=" + to_string(state.n) +
                        " vs dataset n=" + to_string(n) +
                        ". Use the same filters (--keep/--remove/--pheno) in both runs.");
        if (state.x_c != x_c)
            LOGGER.e(0, "Covariate count mismatch: REML state x_c=" + to_string(state.x_c) +
                        " vs current x_c=" + to_string(x_c) +
                        ". Use the same --qcovar/--covar as during --save-reml.");

        // ---- y_adj = y - X*b  (pre-adjusted phenotype) ----
        Eigen::VectorXf y_vec(n);
        for (int i = 0; i < n; ++i) y_vec[i] = static_cast<float>(phenos_vec[i]);
        {
            const Eigen::MatrixXf Xf = X_design.cast<float>();
            y_vec -= Xf * state.b;   // b always loaded because !no_adj_covar
        }

        // ---- Pre-scan setup ----
        Eigen::VectorXf Vi_y(n);
        Eigen::LLT<Eigen::MatrixXf> Vi_llt;

        if (state.is_woodbury) {
            Vi_y = woodbury_apply_Vi_f(state.wb, y_vec);
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
        const bool use_wb = state.is_woodbury;
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

        auto callback = [&](uintptr_t* buf, const vector<uint32_t>& exIdx) {
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
