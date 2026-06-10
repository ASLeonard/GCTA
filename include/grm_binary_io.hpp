#pragma once

#include <Eigen/Dense>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

namespace gcta_grm_io {

// Read a GCTA GRM binary file set (prefix.grm.id, prefix.grm.bin, prefix.grm.N.bin).
// Returns:
//   ids      — "FID\tIID" strings in GRM file order
//   G        — full symmetric n×n matrix (double precision)
//   m_snps   — SNP count from element (0,0) of .grm.N.bin; 0 if N file missing
//
// Implementation notes:
//   - .grm.bin is mmap'd with MADV_SEQUENTIAL for OS read-ahead; avoids a
//     userspace heap allocation of up to several GB at large n.
//   - float→double conversion is batched into a single std::transform pass,
//     keeping buf reads sequential and compiler-vectorisable (AVX2 on Zen 4).
//   - Only the lower triangle of G is filled in the scatter loop; the upper
//     triangle is then materialised via selfadjointView<Lower> to avoid
//     cache-hostile scattered writes into non-sequential columns.
inline void read_grm_binary(const std::string& prefix,
                             std::vector<std::string>& ids,
                             Eigen::MatrixXd& G,
                             double& m_snps)
{
    using std::to_string;

    ids = Pheno::read_sublist(prefix + ".grm.id");
    const int n = static_cast<int>(ids.size());
    if (n == 0) LOGGER.e(0, "GRM id file [" + prefix + ".grm.id] is empty.");

    const size_t tri      = static_cast<size_t>(n) * (n + 1) / 2;
    const size_t byte_len = tri * sizeof(float);

    // ------------------------------------------------------------------ //
    // 1. mmap .grm.bin with MADV_SEQUENTIAL — kernel read-ahead handles   //
    //    prefetching; no userspace buffer needed for the raw float data.   //
    // ------------------------------------------------------------------ //
    const std::string bin_path = prefix + ".grm.bin";
    const int fd = ::open(bin_path.c_str(), O_RDONLY);
    if (fd == -1)
        LOGGER.e(0, "cannot open [" + bin_path + "].");

    {
        struct stat st{};
        if (::fstat(fd, &st) != 0 || static_cast<size_t>(st.st_size) < byte_len) {
            ::close(fd);
            LOGGER.e(0, "unexpected size in [" + bin_path + "]. "
                        "Expected " + to_string(tri) + " float32 entries for n=" +
                        to_string(n) + ".");
        }
    }

    void* raw = ::mmap(nullptr, byte_len, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (raw == MAP_FAILED)
        LOGGER.e(0, "mmap failed for [" + bin_path + "].");

    ::madvise(raw, byte_len, MADV_SEQUENTIAL | MADV_WILLNEED);
    const float* fbuf = static_cast<const float*>(raw);

    // ------------------------------------------------------------------ //
    // 2. Batch convert float→double in one sequential, vectorisable pass. //
    //    Compiler will emit AVX2 vcvtps2pd across the contiguous buffer.  //
    // ------------------------------------------------------------------ //
    std::vector<double> dbuf(tri);
    std::transform(fbuf, fbuf + tri, dbuf.begin(),
                   [](float f) noexcept { return static_cast<double>(f); });

    ::munmap(raw, byte_len);

    // ------------------------------------------------------------------ //
    // 3. Fill only the lower triangle of G (column-major Eigen layout).  //
    //    selfadjointView<Lower> then materialises the upper half in one   //
    //    cache-friendly pass — avoids n²/2 scattered column writes.       //
    // ------------------------------------------------------------------ //
    G.resize(n, n);
    {
        size_t idx = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j <= i; ++j, ++idx)
                G(i, j) = dbuf[idx];
    }
    G = G.selfadjointView<Eigen::Lower>();

    // ------------------------------------------------------------------ //
    // 4. Read SNP count (.grm.N.bin) — only element (0,0) is needed.     //
    // ------------------------------------------------------------------ //
    m_snps = 0.0;
    {
        const std::string n_path = prefix + ".grm.N.bin";
        const int nfd = ::open(n_path.c_str(), O_RDONLY);
        if (nfd != -1) {
            float v0 = 0.0f;
            if (::read(nfd, &v0, sizeof(float)) == sizeof(float))
                m_snps = static_cast<double>(v0);
            ::close(nfd);
        }
        // Missing N file is non-fatal; m_snps remains 0.
    }
}

// Build a mapping from a reference ID list to indices in grm_ids.
// Returns kp[i] = index in grm_ids matching ref_ids[i], or -1 if not found.
//
// Uses unordered_map (O(1) average lookup) rather than std::map (O(log n))
// and reserves capacity upfront to avoid rehashing.
inline std::vector<int> match_ids_to_grm(const std::vector<std::string>& ref_ids,
                                          const std::vector<std::string>& grm_ids)
{
    std::unordered_map<std::string, int> grm_map;
    grm_map.reserve(grm_ids.size() * 2); // factor of 2 keeps load factor ≤ 0.5
    for (int i = 0; i < static_cast<int>(grm_ids.size()); ++i)
        grm_map.emplace(grm_ids[i], i);

    std::vector<int> kp(ref_ids.size(), -1);
    for (int i = 0; i < static_cast<int>(ref_ids.size()); ++i) {
        if (const auto it = grm_map.find(ref_ids[i]); it != grm_map.end())
            kp[i] = it->second;
    }
    return kp;
}

} // namespace gcta_grm_io