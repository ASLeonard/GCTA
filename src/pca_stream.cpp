/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * PCA_stream.cpp — V2 PCA (--pca-stream): computes principal components
 * directly from a GRM file, reusing the same V2 machinery as --mlma-stream's
 * inline REML path (grm_binary_io.hpp / chunked_grm_matvec.hpp /
 * symmetric_eigendecomp.hpp) rather than going through gcta::pca() (V1).
 *
 * Why not just wire --reml-svd-chunked into V1's apply lambda: gcta::pca()
 * routes through manipulate_grm() -> read_grm(), which materializes a dense
 * GRM (twice, if --keep/--remove/--grm-cutoff are used — see
 * manipulate_grm's post-filter resubsetting copy) *before* pca()'s
 * rSVD/Lanczos dispatch is ever reached. Chunking the matvec on top of that
 * would add I/O without reducing peak memory at all. This bypasses that
 * path entirely instead.
 *
 * Scope, deliberately narrower than V1:
 *   - --keep / --remove: supported (plain ID-list intersection/exclusion).
 *   - --grm-cutoff (relatedness pruning): NOT supported. It needs to inspect
 *     off-diagonal GRM values to decide who to drop, which isn't a free
 *     ID-only operation the way --keep/--remove are — a chunked version of
 *     that is separate work. Errors out pointing at V1 instead of guessing.
 *   - --pca-approx is required (rSVD or Lanczos) — V2's only reason to exist
 *     is the chunked/streaming path; use --pca (V1) for exact/dense.
 *
 * Output format matches gcta::pca()'s writer exactly (.eigenval: one value
 * per line; .eigenvec: "FID IID pc1 ... pck", space-separated, no header),
 * so existing consumers — including load_pca_warm_start in MLMA_stream.cpp
 * — keep working regardless of which path produced the file.
 *
 * Plugin wiring follows the same manual-registration pattern as MLMA (see
 * MLMA_stream.cpp / the registers+processMains function-pointer vectors in
 * main): no shared base class, just matching static-method signatures.
 * Add PCAStream::registerOption / PCAStream::processMain to those two
 * vectors to enable it — not done here per your note not to touch that yet.
 */

#include "grm_binary_io.hpp"
#include "chunked_grm_matvec.hpp"
#include "PCA_stream.h"
#include "symmetric_eigendecomp.hpp"

#include <Eigen/Dense>
#include <fstream>
#include <unordered_set>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::to_string;
using std::vector;

map<string, string> PCAStream::options;
map<string, double> PCAStream::options_d;
vector<string> PCAStream::processFunctions;

int PCAStream::registerOption(map<string, vector<string>>& options_in)
{
    // Always capture --out (shared flag)
    // Use the pre-processed "out" key (no dashes), which is always set by main.cpp
    if (!options_in["out"].empty())
        options["out"] = options_in["out"][0];

    const bool has_pca_stream = options_in.find("--pca-stream") != options_in.end();
    if (!has_pca_stream)
        return 0;

    // Only reliable if PCAStream::registerOption runs before whatever else
    // might consume --grm-cutoff for an unrelated purpose (e.g. --make-grm
    // --grm-cutoff in the same invocation) — options_in is shared across all
    // registerOption calls, and this can only see keys not yet erased by an
    // earlier one. Same ordering dependency the whole registers[] list
    // already has; nothing new introduced here.
    if (options_in.find("--grm-cutoff") != options_in.end())
        LOGGER.e(0, "--pca-stream does not support --grm-cutoff yet "
                    "(needs off-diagonal GRM values, not just IDs); use --pca (V1) "
                    "for relatedness pruning, or prune separately first.");

    if (options_in.find("--grm") == options_in.end() || options_in["--grm"].empty())
        LOGGER.e(0, "--pca-stream requires --grm <prefix>.");
    options["grm"] = options_in["--grm"][0];
    options_in.erase("--grm");

    if (options_in.find("--pca-approx") == options_in.end() || options_in["--pca-approx"].empty())
        LOGGER.e(0, "--pca-stream requires --pca-approx <rSVD|Lanczos> "
                    "(V2 has no dense/exact mode — use --pca for that).");
    options["pca_approx"] = options_in["--pca-approx"][0];
    options_in.erase("--pca-approx");

    // TODO: confirm this flag name against whatever V1's --pca uses for its
    // PC-count argument, and reuse that name here instead if it differs.
    if (options_in.find("--pca-out-num") != options_in.end() && !options_in["--pca-out-num"].empty()) {
        options_d["pca_out_num"] = std::stod(options_in["--pca-out-num"][0]);
        options_in.erase("--pca-out-num");
    }

    if (options_in.find("--keep") != options_in.end() && !options_in["--keep"].empty()) {
        options["keep"] = options_in["--keep"][0];
        options_in.erase("--keep");
    }
    if (options_in.find("--remove") != options_in.end() && !options_in["--remove"].empty()) {
        options["remove"] = options_in["--remove"][0];
        options_in.erase("--remove");
    }

    // Same names as --reml-svd-chunked / --reml-svd-chunk-size in
    // MLMA_stream.cpp — same meaning ("read the GRM in tiles instead of
    // densely"), reused verbatim rather than introducing a PCA-specific pair.
    if (options_in.find("--reml-svd-chunked") != options_in.end()) {
        if (!options_in["--reml-svd-chunked"].empty())
            LOGGER.w(0, "--reml-svd-chunked takes no argument; ignoring the supplied value.");
        options["svd_chunked"] = "1";
        options_in.erase("--reml-svd-chunked");
    }
    if (options_in.find("--reml-svd-chunk-size") != options_in.end()) {
        const auto& vals = options_in["--reml-svd-chunk-size"];
        if (vals.empty() || vals[0].empty())
            LOGGER.e(0, "--reml-svd-chunk-size requires a row-count argument.");
        options_d["svd_chunk_size"] = std::stod(vals[0]);
        options_in.erase("--reml-svd-chunk-size");
    }

    processFunctions.push_back("PCAStream");
    options_in.erase("--pca-stream");
    return 1;
}

void PCAStream::processMain()
{
    for (const auto& pf : processFunctions) {
        if (pf != "PCAStream") continue;

        const string grm_pfx    = options.at("grm");
        const string out_prefix = options.at("out");
        const string pca_approx = options.at("pca_approx");

        // ---- Sample set: full GRM ID list, optionally filtered by --keep/--remove ----
        const vector<string> grm_ids = Pheno::read_sublist(grm_pfx + ".grm.id");

        vector<string> analysis_ids = grm_ids;  // preserves GRM file order, same as V1's _keep convention
        if (options.count("keep")) {
            const vector<string> keep_ids = Pheno::read_sublist(options.at("keep"));
            const std::unordered_set<string> keep_set(keep_ids.begin(), keep_ids.end());
            vector<string> filtered;
            filtered.reserve(analysis_ids.size());
            for (const auto& id : analysis_ids)
                if (keep_set.count(id)) filtered.push_back(id);
            analysis_ids = std::move(filtered);
        }
        if (options.count("remove")) {
            const vector<string> remove_ids = Pheno::read_sublist(options.at("remove"));
            const std::unordered_set<string> remove_set(remove_ids.begin(), remove_ids.end());
            vector<string> filtered;
            filtered.reserve(analysis_ids.size());
            for (const auto& id : analysis_ids)
                if (!remove_set.count(id)) filtered.push_back(id);
            analysis_ids = std::move(filtered);
        }
        const int n = static_cast<int>(analysis_ids.size());
        if (n == 0)
            LOGGER.e(0, "--pca-stream: no individuals remain after --keep/--remove filtering.");

        int out_pc_num = options_d.count("pca_out_num")
            ? static_cast<int>(options_d.at("pca_out_num")) : 20;  // TODO: confirm against V1's actual default
        if (out_pc_num <= 0 || out_pc_num > n) out_pc_num = n;
        if (out_pc_num == n)
            LOGGER.e(0, "--pca-stream: requested PC count equals sample size — "
                        "rSVD/Lanczos need out_pc_num < n; use --pca (V1) for a full decomposition.");

        const bool svd_chunked = options.count("svd_chunked") > 0;
        const int  chunk_size  = options_d.count("svd_chunk_size")
            ? static_cast<int>(options_d.at("svd_chunk_size")) : 8000;

        // ---- GRM access: chunked tile reader, or dense (fallback / comparison) ----
        gcta_chunked::TileReader chunked_reader;
        Eigen::MatrixXd G_dense;  // left empty when svd_chunked

        if (svd_chunked) {
            gcta_grm_io::ChunkedGrmHandle handle = gcta_grm_io::make_chunked_grm_reader(grm_pfx, analysis_ids);
            chunked_reader = std::move(handle.reader);
            LOGGER.i(0, "--pca-stream: GRM will be read in " + to_string(chunk_size) +
                        "-row tiles from [" + grm_pfx + "] as needed, not loaded densely.");
        } else {
            vector<string> loaded_ids;
            double m_snps_unused = 0.0;
            Eigen::MatrixXd G_full;
            gcta_grm_io::read_grm_binary(grm_pfx, loaded_ids, G_full, m_snps_unused);
            const vector<int> kp = gcta_grm_io::match_ids_to_grm(analysis_ids, loaded_ids);
            for (int i = 0; i < n; ++i)
                if (kp[i] < 0)
                    LOGGER.e(0, "--pca-stream: individual [" + analysis_ids[i] + "] not found in GRM.");

            // Fast path: no --keep/--remove (or a --keep/--remove that
            // happens to preserve GRM order) means kp is the identity — skip
            // the copy entirely rather than holding G_full and G_dense both
            // fully resident at once for no reason.
            const int n_grm = static_cast<int>(loaded_ids.size());
            bool is_identity = (n_grm == n);
            for (int i = 0; i < n && is_identity; ++i)
                if (kp[i] != i) is_identity = false;

            if (is_identity) {
                G_dense = std::move(G_full);
            } else {
                G_dense.resize(n, n);
                // Column-first traversal for column-major Eigen storage —
                // same pattern as MLMA_stream.cpp's G_n subsetting: for fixed
                // j, varying i reads down one column of G_full (contiguous
                // range, even though kp[i] visits it out of order), rather
                // than jumping across the whole matrix.
                for (int j = 0; j < n; ++j) {
                    const int src_col = kp[j];
                    for (int i = 0; i < n; ++i)
                        G_dense(i, j) = G_full(kp[i], src_col);
                }
                // G_full (n_grm x n_grm, potentially much larger than G_dense
                // when heavily filtered) goes out of scope at the end of this
                // else-block regardless, but free it explicitly right here —
                // no reason to keep it alive through the rest of this branch.
                G_full.resize(0, 0);
            }
        }

        auto apply = [&](const auto& X) -> Eigen::MatrixXd {
            if (svd_chunked)
                return gcta_chunked::chunked_symmetric_matvec(chunked_reader, n, chunk_size, X);
            return G_dense * X;
        };

        gcta_eigh::EighResult res;
        try {
            if (pca_approx == "rSVD") {
                const int oversample = gcta_eigh::recommended_oversample(out_pc_num);
                res = gcta_eigh::randomized_symmetric_eigh(apply, n, out_pc_num, oversample, 3);
            } else if (pca_approx == "Lanczos") {
                const int ncv = std::min(n, std::max(3 * out_pc_num + 1, 30));
                res = gcta_eigh::lanczos_symmetric_eigh(apply, n, out_pc_num, ncv);
            } else {
                LOGGER.e(0, "--pca-approx: unrecognised method '" + pca_approx + "'. Use 'rSVD' or 'Lanczos'.");
            }
        } catch (const std::exception& e) {
            LOGGER.e(0, string("--pca-stream failed: ") + e.what());
        }

        // Neither is needed past this point — G_dense (up to n x n) and the
        // chunked reader's mmap can both be released before the O(n) file
        // writes below, rather than sitting at peak size through them.
        G_dense.resize(0, 0);
        chunked_reader = nullptr;

        // ---- Write .eigenval / .eigenvec — identical format to gcta::pca()'s writer ----
        const string eval_file = out_prefix + ".eigenval";
        std::ofstream o_eval(eval_file);
        if (!o_eval) LOGGER.e(0, "cannot open the file [" + eval_file + "] to write.");
        for (int i = 0; i < out_pc_num; ++i) o_eval << res.eigenvalues(i) << "\n";
        o_eval.close();
        LOGGER.i(0, "Eigenvalues of " + to_string(n) + " individuals have been saved in [" + eval_file + "].");

        const string evec_file = out_prefix + ".eigenvec";
        std::ofstream o_evec(evec_file);
        if (!o_evec) LOGGER.e(0, "cannot open the file [" + evec_file + "] to write.");

        // res.eigenvectors is n x out_pc_num, column-major: reading it row by
        // row (one individual's full PC vector at a time, which is what the
        // output format needs) strides by n elements between each PC value.
        // Materializing the transpose once (out_pc_num x n, still
        // column-major) makes each individual's PCs a contiguous column
        // instead — one O(n*out_pc_num) sequential copy up front in exchange
        // for a cache-friendly access pattern in the n-times-larger write
        // loop that follows.
        const Eigen::MatrixXd evec_t = res.eigenvectors.transpose();
        for (int i = 0; i < n; ++i) {
            // analysis_ids[i] is "FID\tIID" (Pheno::read_sublist convention,
            // same as match_ids_to_grm's expected input elsewhere); V1's
            // writer joins FID/IID with a space, so convert to match exactly.
            string id = analysis_ids[i];
            const size_t tab = id.find('\t');
            if (tab != string::npos) id[tab] = ' ';
            o_evec << id;
            for (int j = 0; j < out_pc_num; ++j) o_evec << " " << evec_t(j, i);
            o_evec << "\n";
        }
        o_evec.close();
        LOGGER.i(0, to_string(n) + " individuals, " + to_string(out_pc_num) +
                    " eigenvectors saved in [" + evec_file + "].");
    }
}
