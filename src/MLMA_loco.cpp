/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * MLMALoco — v2 LOCO MLMA.  Per-chromosome REML via reml::compute() + streaming
 * association test via Geno::loopDouble.
 *
 * Mandatory flags:
 *   --mlma-loco-stream    (trigger)
 *   --loco-manifest <f>   (TSV: chrom, bfile_prefix, grm_chr_prefix)
 *   --grm <prefix>        (all-autosome GRM)
 *   --pheno <f>
 *
 * Optional:
 *   --qcovar / --covar
 *   --reml-woodbury [k]
 *   --reml-trace-approx [n]
 *   --log-pval
 *   --out <prefix>
 */

#include "MLMA_loco.h"
#include "RemlCtx.hpp"
#include "RemlEngine.hpp"
#include "RemlState.hpp"
#include "mlma_woodbury.hpp"
#include "Logger.h"
#include "cpu.h"
#include "Pheno.h"
#include "Marker.h"
#include "Geno.h"
#include "Covar.h"

#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>

using std::map;
using std::string;
using std::vector;
using std::to_string;
using std::function;

// ---------------------------------------------------------------------------
// Static member definitions
// ---------------------------------------------------------------------------
map<string, string>  MLMALoco::options;
map<string, double>  MLMALoco::options_d;
vector<string>       MLMALoco::processFunctions;

// ---------------------------------------------------------------------------
// GRM binary I/O helpers (file-local)
// ---------------------------------------------------------------------------
namespace {

// Read a GCTA GRM binary file set (prefix.grm.id, prefix.grm.bin, prefix.grm.N.bin).
// Returns:
//   ids      — "FID\tIID" strings in GRM file order
//   G        — full symmetric n×n matrix (double precision)
//   m_snps   — SNP count from diagonal (0,0) of .grm.N.bin; 0 if N file missing
void read_grm_binary(const string& prefix,
                     vector<string>& ids,
                     Eigen::MatrixXd& G,
                     double& m_snps)
{
    // Read IDs
    ids = Pheno::read_sublist(prefix + ".grm.id");
    const int n = static_cast<int>(ids.size());
    if (n == 0) LOGGER.e(0, "GRM id file [" + prefix + ".grm.id] is empty.");

    const size_t tri = static_cast<size_t>(n) * (n + 1) / 2;

    // Read GRM matrix (.grm.bin — float32, lower triangle row-by-row)
    {
        std::ifstream f(prefix + ".grm.bin", std::ios::binary);
        if (!f.is_open()) LOGGER.e(0, "cannot open [" + prefix + ".grm.bin].");
        vector<float> buf(tri);
        f.read(reinterpret_cast<char*>(buf.data()),
               static_cast<std::streamsize>(tri * sizeof(float)));
        if (static_cast<size_t>(f.gcount()) != tri * sizeof(float))
            LOGGER.e(0, "unexpected EOF in [" + prefix + ".grm.bin]. "
                        "Expected " + to_string(tri) + " float32 entries for n=" + to_string(n) + ".");
        G.resize(n, n);
        size_t idx = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j <= i; ++j, ++idx) {
                G(i, j) = static_cast<double>(buf[idx]);
                G(j, i) = static_cast<double>(buf[idx]);
            }
    }

    // Read SNP count matrix (.grm.N.bin — float32, same layout)
    m_snps = 0.0;
    {
        std::ifstream fn(prefix + ".grm.N.bin", std::ios::binary);
        if (fn.is_open()) {
            float v0 = 0.0f;
            fn.read(reinterpret_cast<char*>(&v0), sizeof(float));
            if (fn.gcount() == sizeof(float)) m_snps = static_cast<double>(v0);
        }
    }
}

// Build a mapping from a reference ID list (ref_ids) to indices in grm_ids.
// Returns kp[i] = index in grm_ids matching ref_ids[i], or -1 if not found.
vector<int> match_ids_to_grm(const vector<string>& ref_ids,
                              const vector<string>& grm_ids)
{
    map<string, int> grm_map;
    for (int i = 0; i < (int)grm_ids.size(); ++i)
        grm_map[grm_ids[i]] = i;
    vector<int> kp(ref_ids.size(), -1);
    for (int i = 0; i < (int)ref_ids.size(); ++i) {
        auto it = grm_map.find(ref_ids[i]);
        if (it != grm_map.end()) kp[i] = it->second;
    }
    return kp;
}

// Manifest row
struct ManifestRow {
    int    chrom;
    string bfile;
    string grm_chr;
};

vector<ManifestRow> read_manifest(const string& filename)
{
    std::ifstream f(filename);
    if (!f.is_open()) LOGGER.e(0, "cannot open manifest file [" + filename + "].");
    vector<ManifestRow> rows;
    string line;
    int lineno = 0;
    while (std::getline(f, line)) {
        ++lineno;
        // strip \r
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        ManifestRow r;
        string chrom_str;
        if (!(ss >> chrom_str >> r.bfile >> r.grm_chr)) {
            LOGGER.e(0, "manifest [" + filename + "] line " + to_string(lineno)
                      + ": expected 3 columns (chrom, bfile, grm_chr).");
        }
        r.chrom = std::stoi(chrom_str);
        rows.push_back(std::move(r));
    }
    if (rows.empty()) LOGGER.e(0, "manifest file [" + filename + "] has no data rows.");
    return rows;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// registerOption
// ---------------------------------------------------------------------------
int MLMALoco::registerOption(map<string, vector<string>>& options_in)
{
    // Common option: --out
    if (options_in.find("out") != options_in.end() && !options_in["out"].empty())
        options["out"] = options_in["out"][0];

    const bool has_flag = (options_in.find("--mlma-loco-stream") != options_in.end());
    if (!has_flag) return 0;

    // --loco-manifest (mandatory)
    if (options_in.find("--loco-manifest") == options_in.end()
            || options_in["--loco-manifest"].empty())
        LOGGER.e(0, "--mlma-loco-stream requires --loco-manifest <file>.");
    options["manifest"] = options_in["--loco-manifest"][0];
    options_in.erase("--loco-manifest");

    // --grm (mandatory — all-chromosome GRM)
    if (options_in.find("--grm") != options_in.end() && !options_in["--grm"].empty()) {
        options["grm_all"] = options_in["--grm"][0];
        options_in.erase("--grm");
    } else {
        LOGGER.e(0, "--mlma-loco-stream requires --grm <prefix> (all-autosome GRM).");
    }

    // Optional flags
    if (options_in.find("--log-pval") != options_in.end()) {
        options["log_pval"] = "1";
        options_in.erase("--log-pval");
    }

    // --reml-woodbury [k]  (k optional; -1 = auto-k)
    if (options_in.find("--reml-woodbury") != options_in.end()) {
        const auto& vals = options_in["--reml-woodbury"];
        if (!vals.empty() && !vals[0].empty())
            options_d["woodbury_rank"] = std::stod(vals[0]);
        else
            options_d["woodbury_rank"] = -1.0;   // auto-k
        options_in.erase("--reml-woodbury");
    }

    // --reml-trace-approx [n]
    if (options_in.find("--reml-trace-approx") != options_in.end()) {
        const auto& vals = options_in["--reml-trace-approx"];
        options["trace_approx"] = "1";
        if (!vals.empty() && !vals[0].empty())
            options_d["trace_approx_nprobes"] = std::stod(vals[0]);
        options_in.erase("--reml-trace-approx");
    }

    // --reml-maxit
    if (options_in.find("--reml-maxit") != options_in.end()
            && !options_in["--reml-maxit"].empty()) {
        options_d["reml_maxit"] = std::stod(options_in["--reml-maxit"][0]);
        options_in.erase("--reml-maxit");
    }

    processFunctions.push_back("MLMALoco");
    options_in.erase("--mlma-loco-stream");
    return 1;
}

// ---------------------------------------------------------------------------
// processMain
// ---------------------------------------------------------------------------
void MLMALoco::processMain()
{
    for (const auto& pf : processFunctions) {
        if (pf != "MLMALoco") continue;

        const string manifest_file = options.at("manifest");
        const string grm_all_pfx  = options.at("grm_all");
        const string out_prefix   = options.at("out");
        const bool   log_pval     = options.count("log_pval") > 0;

        // REML config
        const int    woodbury_rank = (options_d.count("woodbury_rank") > 0)
            ? static_cast<int>(options_d.at("woodbury_rank")) : 0;
        const bool   trace_approx = (options.count("trace_approx") > 0);
        const int    trace_nprobes = (options_d.count("trace_approx_nprobes") > 0)
            ? static_cast<int>(options_d.at("trace_approx_nprobes")) : 90;
        const int    reml_maxit = (options_d.count("reml_maxit") > 0)
            ? static_cast<int>(options_d.at("reml_maxit")) : 100;

        // ---- 1. Read manifest ----
        vector<ManifestRow> manifest = read_manifest(manifest_file);
        LOGGER << "[LOCO] Manifest: " << manifest.size() << " chromosome(s)." << std::endl;

        // ---- 2. Load all-chr GRM ----
        LOGGER << "\n[LOCO] Loading all-autosome GRM from [" << grm_all_pfx << "] ..." << std::endl;
        vector<string> grm_all_ids;
        Eigen::MatrixXd G_all;
        double m_all = 0.0;
        read_grm_binary(grm_all_pfx, grm_all_ids, G_all, m_all);
        if (m_all <= 0.0)
            LOGGER.e(0, "[LOCO] Could not read marker count from [" + grm_all_pfx
                      + ".grm.N.bin]. Rebuild GRM without --no-grm-N.");
        LOGGER << "[LOCO] G_all: " << grm_all_ids.size() << " individuals, "
               << static_cast<long long>(m_all) << " markers." << std::endl;

        // ---- 3. Load phenotype ----
        Pheno* pheno  = new Pheno();
        const uint32_t n_pheno_init = pheno->count_keep();
        vector<string> pheno_ids;
        vector<double> phenos_vec;
        pheno->get_pheno(pheno_ids, phenos_vec);
        if (pheno_ids.empty())
            LOGGER.e(0, "--mlma-loco-stream requires --pheno to be specified.");

        // ---- 4. Load covariates ----
        Covar covar;
        bool  has_covar = false;
        int   x_c       = 1;    // intercept only by default
        Eigen::MatrixXd X_design_full;  // n_pheno × x_c (before GRM intersection)

        {
            vector<uint32_t> remain_index, covar_index;
            has_covar = covar.getCommonSampleIndex(pheno_ids, remain_index, covar_index);

            if (has_covar) {
                LOGGER.i(0, to_string(remain_index.size()) +
                            " overlapping individuals with non-missing covariate data.");
                pheno->filter_keep_index(remain_index);
                // reorder phenos_vec
                vector<double> phenos2(remain_index.size());
                for (size_t i = 0; i < remain_index.size(); ++i)
                    phenos2[i] = phenos_vec[remain_index[i]];
                phenos_vec = std::move(phenos2);
            } else if (covar.hasCovar()) {
                LOGGER.e(0, "no overlapping individuals with non-missing covariate data.");
            } else {
                remain_index.resize(pheno_ids.size());
                std::iota(remain_index.begin(), remain_index.end(), 0u);
            }

            const int n_new = static_cast<int>(pheno->count_keep());
            // Rebuild pheno_ids after filtering
            pheno_ids.clear();
            phenos_vec.resize(n_new);
            {
                vector<string> all_ids_new;
                vector<double> all_pheno_new;
                pheno->get_pheno(all_ids_new, all_pheno_new);
                pheno_ids = std::move(all_ids_new);
                phenos_vec = std::move(all_pheno_new);
            }

            if (has_covar) {
                vector<string> final_ids = pheno->get_id(0, n_new - 1, "\t");
                vector<double> covar_X;
                vector<uint32_t> ck_dummy;
                covar.getCovarX(final_ids, covar_X, ck_dummy);
                const int covar_cols = (n_new > 0)
                    ? static_cast<int>(covar_X.size()) / n_new : 0;
                x_c = 1 + covar_cols;
                X_design_full.resize(n_new, x_c);
                X_design_full.col(0).setOnes();
                if (covar_cols > 0)
                    X_design_full.rightCols(covar_cols) =
                        Eigen::Map<const Eigen::MatrixXd>(covar_X.data(), n_new, covar_cols);
            } else {
                const int n_new2 = static_cast<int>(pheno->count_keep());
                X_design_full.resize(n_new2, 1);
                X_design_full.setOnes();
            }
        }

        // ---- 5. Intersect pheno+covar individuals with GRM ----
        // Match pheno IDs (FID\tIID) to GRM IDs (FID\tIID from read_sublist)
        vector<int> kp = match_ids_to_grm(pheno_ids, grm_all_ids);
        // Build the final keep list: only individuals present in GRM
        vector<int> final_keep;
        final_keep.reserve(pheno_ids.size());
        for (int i = 0; i < (int)pheno_ids.size(); ++i)
            if (kp[i] >= 0) final_keep.push_back(i);

        if (final_keep.empty())
            LOGGER.e(0, "[LOCO] No individuals overlap between phenotype and GRM. "
                        "Check that GRM was built from the same sample.");

        const int n = static_cast<int>(final_keep.size());
        LOGGER << "[LOCO] After GRM intersection: " << n << " individuals for analysis." << std::endl;

        // Subset G_all to the n analysis individuals.
        // Fast path: if GRM sample == analysis sample (in order), move instead of copy.
        Eigen::MatrixXd G_all_n;
        {
            const int n_grm = static_cast<int>(grm_all_ids.size());
            vector<int> grm_rows(n);
            bool is_identity = (n_grm == n);
            for (int i = 0; i < n; ++i) {
                grm_rows[i] = kp[final_keep[i]];
                if (grm_rows[i] != i) is_identity = false;
            }
            if (is_identity) {
                G_all_n = std::move(G_all);  // O(1): GRM sample == analysis sample
            } else {
                G_all_n.resize(n, n);
                // Column-first traversal: accesses contiguous columns of G_all
                // (column-major Eigen storage), reducing cache misses by ~n×.
                for (int j = 0; j < n; ++j) {
                    const int src_col = grm_rows[j];
                    for (int i = 0; i < n; ++i)
                        G_all_n(i, j) = G_all(grm_rows[i], src_col);
                }
                G_all.resize(0, 0);
            }
        }

        // Subset y and X
        Eigen::VectorXd y_vec(n);
        Eigen::MatrixXd X_mat(n, x_c);
        for (int i = 0; i < n; ++i) {
            y_vec[i] = phenos_vec[final_keep[i]];
            X_mat.row(i) = X_design_full.row(final_keep[i]);
        }
        X_design_full.resize(0, 0);  // free

        // Also build the final analysis IDs (for GRM subset matching per chr)
        vector<string> analysis_ids(n);
        for (int i = 0; i < n; ++i) analysis_ids[i] = pheno_ids[final_keep[i]];

        // Filter pheno to exactly these individuals (needed by Geno::loopDouble)
        {
            // Rebuild keep indices inside Pheno to match analysis_ids order
            // The simplest way: get the raw keep indices from Pheno and re-filter
            // We create a new Pheno object per chromosome bfile anyway, so this
            // Pheno object is only needed for reference — no further filtering needed.
        }

        // ---- 6. Open output file ----
        const string out_file = out_prefix + ".loco.mlma";
        LOGGER << "\n[LOCO] LOCO results will be streamed to [" << out_file << "] ..." << std::endl;
        std::ofstream ofile(out_file);
        if (!ofile) LOGGER.e(0, "cannot open [" + out_file + "] for writing.");
        ofile << "Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\t"
              << (log_pval ? "log_p" : "p") << "\n";

        // ---- 7. Per-chromosome LOCO loop ----
        LOGGER << "\n[LOCO] Running LOCO MLMA over " << manifest.size()
               << " chromosome(s) from manifest." << std::endl;
        if (woodbury_rank != 0)
            LOGGER << "[LOCO] --reml-woodbury active; basis warm-started between chromosomes." << std::endl;
        else
            LOGGER << "[LOCO] Tip: add --reml-woodbury <k> for faster per-chr REML." << std::endl;

        // Woodbury basis for warm-starting across chromosomes
        Eigen::MatrixXd Uk_warmstart;
        Eigen::VectorXd dk_warmstart;

        // I/O batch buffer to reduce write() syscall overhead
        std::vector<char> io_buf(4 << 20);  // 4 MiB
        size_t io_pos = 0;
        auto flush_io = [&]() {
            if (io_pos > 0) {
                ofile.write(io_buf.data(), static_cast<std::streamsize>(io_pos));
                io_pos = 0;
            }
        };

        for (size_t ci = 0; ci < manifest.size(); ++ci) {
            const ManifestRow& row = manifest[ci];
            LOGGER << "\n-----------------------------------\n#Chr " << row.chrom << ":" << std::endl;

            // ---- 7a. Load G_chr ----
            vector<string> grm_chr_ids;
            Eigen::MatrixXd G_chr_raw;
            double m_chr = 0.0;
            LOGGER << "  Reading per-chr GRM from [" << row.grm_chr << "] ..." << std::endl;
            read_grm_binary(row.grm_chr, grm_chr_ids, G_chr_raw, m_chr);
            if (m_chr <= 0.0)
                LOGGER.e(0, "[LOCO] Chr " + to_string(row.chrom)
                          + ": could not read marker count from [" + row.grm_chr
                          + ".grm.N.bin]. Rebuild per-chr GRM without --no-grm-N.");

            // Match analysis_ids to chr GRM IDs
            vector<int> kp_chr = match_ids_to_grm(analysis_ids, grm_chr_ids);
            for (int i = 0; i < n; ++i)
                if (kp_chr[i] < 0)
                    LOGGER.e(0, "[LOCO] Chr " + to_string(row.chrom)
                              + ": individual [" + analysis_ids[i]
                              + "] not found in per-chr GRM [" + row.grm_chr + "].");

            // Subset G_chr to analysis individuals.
            // Fast path: if chr-GRM sample == analysis sample (in order), move.
            Eigen::MatrixXd G_chr;
            {
                const int n_grm_chr = static_cast<int>(grm_chr_ids.size());
                bool chr_is_identity = (n_grm_chr == n);
                for (int i = 0; i < n && chr_is_identity; ++i)
                    if (kp_chr[i] != i) chr_is_identity = false;

                if (chr_is_identity) {
                    G_chr = std::move(G_chr_raw);
                } else {
                    G_chr.resize(n, n);
                    // Column-first traversal for column-major Eigen storage.
                    for (int j = 0; j < n; ++j) {
                        const int src_col = kp_chr[j];
                        for (int i = 0; i < n; ++i)
                            G_chr(i, j) = G_chr_raw(kp_chr[i], src_col);
                    }
                    G_chr_raw.resize(0, 0);
                }
            }

            // ---- 7b. Compute LOCO GRM in-place via G_chr ----
            // G_loco = (G_all_n * m_all - G_chr * m_chr) / m_loco
            // Reuse G_chr's buffer to avoid a 3rd n×n allocation (saves ~1.8 GB peak).
            const double m_loco = m_all - m_chr;
            if (m_loco <= 0.0)
                LOGGER.e(0, "[LOCO] Chr " + to_string(row.chrom) + ": m_chr=" +
                            to_string(static_cast<long long>(m_chr))
                          + " >= m_all=" + to_string(static_cast<long long>(m_all)) + ".");
            G_chr.array() *= (-m_chr / m_loco);
            G_chr.noalias() += (m_all / m_loco) * G_all_n;
            // G_chr now holds G_loco; G_chr_raw is already freed.

            // ---- 7c. Fill RemlCtx ----
            RemlCtx ctx;
            ctx.n   = n;
            ctx.X_c = x_c;
            ctx.X   = X_mat;
            ctx.y   = y_vec;
            ctx.A.resize(2);
            ctx.A[0] = std::move(G_chr);  // consumed by compute_woodbury_basis or used as-is
            // ctx.A[1] left default (size 0 == identity convention for residual)
            ctx.r_indx = {0, 1};
            ctx.var_name = {"V(G)", "V(e)"};
            ctx.hsq_name = {"V(G)/Vp"};
            ctx.out = out_prefix + ".chr" + to_string(row.chrom);

            // REML config
            ctx.reml_mtd       = 0;  // AI-REML
            ctx.reml_max_iter  = reml_maxit;
            ctx.reml_inv_mtd   = 0;  // LLT
            ctx.woodbury_rank  = woodbury_rank;
            ctx.woodbury_buffer_factor = 1.5;
            ctx.reml_trace_approx = trace_approx;
            ctx.reml_trace_approx_nprobes = trace_nprobes;

            // Warm-start woodbury basis from previous chr
            if (woodbury_rank != 0 && Uk_warmstart.rows() > 0) {
                ctx.Uk = Uk_warmstart;
                ctx.dk = dk_warmstart;
            }

            // ---- 7d. Run REML ----
            reml::compute(ctx, /*priors=*/{}, /*priors_var=*/{}, /*no_constrain=*/false);

            // Save eigenvectors for warm-start in next chromosome
            if (ctx.Vi_use_woodbury && ctx.Uk.rows() > 0) {
                Uk_warmstart = ctx.Uk;
                dk_warmstart = ctx.dk;
            }

            LOGGER << "[LOCO] Chr " << row.chrom << " variance components: ";
            for (double vc : ctx.varcmp) LOGGER << vc << " ";
            LOGGER << std::endl;

            // ---- 7e. Build RemlState ----
            RemlState rs = reml::build_reml_state(ctx);

            // ---- 7f. Pre-adjusted phenotype: y_adj = y - X*b ----
            Eigen::VectorXf y_adj(n);
            {
                Eigen::VectorXd b_d = ctx.b;
                Eigen::VectorXd y_adj_d = y_vec - X_mat * b_d;
                for (int i = 0; i < n; ++i) y_adj[i] = static_cast<float>(y_adj_d[i]);
            }

            // ---- 7g. Compute Vi_y_adj ----
            Eigen::VectorXf Vi_y(n);
            Eigen::LLT<Eigen::MatrixXf> Vi_llt;
            const bool use_wb  = rs.is_woodbury;
            const bool use_llt = rs.is_llt;

            if (use_wb) {
                Vi_y = woodbury_apply_Vi_f(rs.wb, y_adj);
            } else if (use_llt) {
                // V^{-1} y = L^{-T}(L^{-1} y) via two float triangular solves.
                // Avoids dpotri + second Cholesky (saves 2×O(n³/3) per chromosome).
                Vi_y = y_adj;
                cblas_strsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                            n, rs.Vi_L_f.data(), n, Vi_y.data(), 1);
                cblas_strsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit,
                            n, rs.Vi_L_f.data(), n, Vi_y.data(), 1);
            } else {
                Vi_y.noalias() = rs.Vi * y_adj;
                Vi_llt = Eigen::LLT<Eigen::MatrixXf>(rs.Vi);
                if (Vi_llt.info() != Eigen::Success)
                    LOGGER.e(0, "[LOCO] Chr " + to_string(row.chrom)
                              + ": Vi matrix is not positive definite.");
                rs.Vi.resize(0, 0);
            }

            // ---- 7h. Stream SNPs via Geno::loopDouble ----
            Pheno*  pheno_chr  = new Pheno();
            Marker* marker_chr = new Marker();
            Geno*   geno_chr   = new Geno(pheno_chr, marker_chr);

            const uint32_t total_m = marker_chr->count_extract();
            LOGGER << "[LOCO] Chr " << row.chrom << ": " << total_m << " SNPs to test." << std::endl;

            constexpr int BLOCK = 10000;
            Eigen::MatrixXf X_block(n, BLOCK);
            Eigen::VectorXf Xt_Vi_y(BLOCK);
            Eigen::VectorXf xvx_diag(BLOCK);
            vector<GenoBufItem> gbuf_items(BLOCK);
            vector<uint8_t>     valid_v(BLOCK, 0);
            vector<float>       af_v(BLOCK, 0.0f);

            int snp_done = 0;
            int last_pct = -1;
            const WoodburyMLMACache& wb = rs.wb;

            auto callback = [&](uintptr_t* buf, const vector<uint32_t>& exIdx) {
                const int bs = static_cast<int>(exIdx.size());

                for (int i = 0; i < bs; ++i) {
                    valid_v[i] = 0;
                    af_v[i]    = 0.0f;
                    GenoBufItem& item = gbuf_items[i];
                    item.extractedMarkerIndex = exIdx[i];
                    geno_chr->getGenoDouble(buf, i, &item);
                    if (!item.valid) { X_block.col(i).setZero(); continue; }
                    valid_v[i] = 1;
                    af_v[i]    = static_cast<float>(item.af);
                    X_block.col(i) =
                        Eigen::Map<const Eigen::VectorXd>(item.geno.data(), n).cast<float>();
                }

                Xt_Vi_y.head(bs).noalias() =
                    X_block.leftCols(bs).transpose() * Vi_y;

                if (use_wb) {
                    woodbury_xvx_diag_block(wb, X_block, bs, xvx_diag);
                } else if (use_llt) {
                    // x^T V^{-1} x = ||L^{-1} x||^2 (L = lower Cholesky of V).
                    // In-place STRSM overwrites X_block with L^{-1} X_block.
                    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                                n, bs, 1.0f,
                                rs.Vi_L_f.data(), n,
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

                for (int i = 0; i < bs; ++i) {
                    const uint32_t raw = marker_chr->getRawIndex(exIdx[i]);
                    const std::string& chr_s = marker_chr->getRawChr(raw);
                    const std::string& name  = marker_chr->getRawName(raw);
                    const unsigned     bp    = static_cast<unsigned>(marker_chr->getRawBp(raw));
                    const std::string& a1    = marker_chr->getRawA1(raw);
                    const std::string& a2    = marker_chr->getRawA2(raw);

                    char line[512];
                    int len;
                    float beta_val, se_val;
                    double pval_val;
                    if (!valid_v[i] || !mlma_snp_stat(Xt_Vi_y[i], xvx_diag[i], log_pval,
                                                      beta_val, se_val, pval_val)) {
                        len = std::snprintf(line, sizeof(line),
                                            "%s\t%s\t%u\t%s\t%s\tNA\tNA\tNA\tNA\n",
                                            chr_s.c_str(), name.c_str(), bp,
                                            a1.c_str(), a2.c_str());
                    } else {
                        len = std::snprintf(line, sizeof(line),
                                            "%s\t%s\t%u\t%s\t%s\t%.6g\t%.6g\t%.6g\t%.6g\n",
                                            chr_s.c_str(), name.c_str(), bp,
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

            const vector<uint32_t>& extractIndex = marker_chr->get_extract_index();
            geno_chr->loopDouble(extractIndex, BLOCK,
                                 /*bMakeGeno*/   true,
                                 /*bGenoCenter*/ true,
                                 /*bGenoStd*/    false,
                                 /*bMakeMiss*/   true,
                                 {callback});

            flush_io();
            LOGGER << "[LOCO] Chr " << row.chrom << " done." << std::endl;

            delete geno_chr;
            delete marker_chr;
            delete pheno_chr;
        }

        flush_io();
        LOGGER << "\n[LOCO] Results saved to [" << out_file << "]." << std::endl;
        ofile.close();

        delete pheno;
    }
}
