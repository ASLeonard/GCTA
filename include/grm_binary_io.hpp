/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * grm_binary_io.hpp — Shared inline helpers for reading GCTA GRM binary files.
 *
 * Declared inline so the definitions can be included in multiple translation
 * units (including unity-build units) without ODR violations.
 *
 * Used by:  src/MLMA_loco.cpp  (per-chromosome + all-chromosome GRM loading)
 *           src/MLMA_stream.cpp (all-chromosome GRM for inline REML)
 */
#pragma once

#include "Logger.h"
#include "Pheno.h"
#include <Eigen/Dense>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace gcta_grm_io {

// Read a GCTA GRM binary file set (prefix.grm.id, prefix.grm.bin, prefix.grm.N.bin).
// Returns:
//   ids      — "FID\tIID" strings in GRM file order
//   G        — full symmetric n×n matrix (double precision)
//   m_snps   — SNP count from element (0,0) of .grm.N.bin; 0 if N file is missing
inline void read_grm_binary(const std::string& prefix,
                             std::vector<std::string>& ids,
                             Eigen::MatrixXd& G,
                             double& m_snps)
{
    using std::to_string;

    ids = Pheno::read_sublist(prefix + ".grm.id");
    const int n = static_cast<int>(ids.size());
    if (n == 0) LOGGER.e(0, "GRM id file [" + prefix + ".grm.id] is empty.");

    const size_t tri = static_cast<size_t>(n) * (n + 1) / 2;

    // Read GRM matrix (.grm.bin — float32, lower triangle row-by-row)
    {
        std::ifstream f(prefix + ".grm.bin", std::ios::binary);
        if (!f.is_open()) LOGGER.e(0, "cannot open [" + prefix + ".grm.bin].");
        std::vector<float> buf(tri);
        f.read(reinterpret_cast<char*>(buf.data()),
               static_cast<std::streamsize>(tri * sizeof(float)));
        if (static_cast<size_t>(f.gcount()) != tri * sizeof(float))
            LOGGER.e(0, "unexpected EOF in [" + prefix + ".grm.bin]. "
                        "Expected " + to_string(tri) + " float32 entries for n=" +
                        to_string(n) + ".");
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

// Build a mapping from a reference ID list to indices in grm_ids.
// Returns kp[i] = index in grm_ids matching ref_ids[i], or -1 if not found.
inline std::vector<int> match_ids_to_grm(const std::vector<std::string>& ref_ids,
                                          const std::vector<std::string>& grm_ids)
{
    std::map<std::string, int> grm_map;
    for (int i = 0; i < static_cast<int>(grm_ids.size()); ++i)
        grm_map[grm_ids[i]] = i;
    std::vector<int> kp(ref_ids.size(), -1);
    for (int i = 0; i < static_cast<int>(ref_ids.size()); ++i) {
        auto it = grm_map.find(ref_ids[i]);
        if (it != grm_map.end()) kp[i] = it->second;
    }
    return kp;
}

} // namespace gcta_grm_io
