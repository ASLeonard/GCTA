#pragma once
#include <map>
#include <string>
#include <vector>

// MLMALoco — v2 LOCO (Leave-One-Chromosome-Out) MLMA analysis.
//
// Usage:
//   gcta64 --mlma-loco-stream --loco-manifest <file> --grm <prefix>
//          --pheno <file> [--qcovar <file>] [--covar <file>]
//          [--reml-woodbury [k]] [--reml-trace-approx [n]]
//          [--log-pval] [--mlma-no-preadj-covar] --out <prefix>
//
// Manifest TSV format (no header):
//   chrom  bfile_prefix  grm_chr_prefix
//
// For each chromosome, the LOCO GRM is computed as:
//   G_loco = (G_all * m_all - G_chr * m_chr) / (m_all - m_chr)
// REML is run on G_loco, then SNPs from bfile are tested.

class MLMALoco {
public:
    static int  registerOption(std::map<std::string, std::vector<std::string>>& options_in);
    static void processMain();

private:
    static std::map<std::string, std::string> options;
    static std::map<std::string, double>      options_d;
    static std::map<std::string, std::vector<double>> options_vd;
    static std::vector<std::string>           processFunctions;
};
