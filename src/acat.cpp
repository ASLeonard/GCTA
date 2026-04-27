#include <cmath>
#include <ctime>
#include <cctype>
#include <iostream>
#include <iomanip>

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#ifndef _WIN32
#include <gsl/gsl_sf_gamma.h>
#else
#include <math.h>
#define gsl_sf_gamma tgamma
#endif

#if defined _WIN32 || defined _WIN64
#include <corecrt_math_defines.h>
#endif

//#define ACAT_UTILITY

// ---------------------------------------------------------------------------
// Macros (math core — unchanged from original)
// ---------------------------------------------------------------------------

#define BETA_DENSITY(x, a, b) \
    (((gsl_sf_gamma(a + b)) / (gsl_sf_gamma(a) * gsl_sf_gamma(b))) \
     * pow((x), a - 1.0) * pow(1.0 - (x), b - 1.0))

#define PCAUCHY(x) \
    ((1.0L / M_PI) * atan((x - 0L) / 1.0L) + 0.5L)

#define PCAUCHY2(x) \
    (atan(1.0 / (x)) / M_PI)

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

struct SnpNode {
    unsigned long snp_pos  = 0;
    double        weight   = 0.0;
    double        p_value  = 0.0;
};

struct GeneNode {
    unsigned long gene_start = 0;
    unsigned long gene_end   = 0;
    std::string   gene_name;
    std::vector<SnpNode *> snp_ptrs;   // replaces SNP_POINTER_NODE linked list
};

// Per-chromosome bucket.
// gene_order_list encodes (start, end, running_max_right) triples — same
// semantics as the original flat array, built after sorting.
struct ChromBucket {
    std::vector<GeneNode>  genes;      // sorted by gene_start after loading
    std::vector<SnpNode>   snps;
    std::vector<unsigned long> gene_order_list;  // triples: [start, end, max_right]
};

// The top-level table is now a plain map from chromosome string → bucket.
using HashTable = std::unordered_map<std::string, ChromBucket>;

// ---------------------------------------------------------------------------
// Chromosome key helpers  (replaces hash_func / unhash_func)
// ---------------------------------------------------------------------------

// Normalise the raw chromosome field to a canonical string key.
// Accepts numeric strings ("1"–"22"), "X", "Y", "XY".
static std::string chrom_key(const std::string & s)
{
    if (s == "X" || s == "Y" || s == "XY")
        return s;
    if (!s.empty() && std::isdigit(static_cast<unsigned char>(s[0])))
        return s;           // keep as-is (leading zeros would differ, but input is clean)
    throw std::runtime_error("Chromosome not recognised: " + s);
}

// ---------------------------------------------------------------------------
// File-parsing helpers
// ---------------------------------------------------------------------------

// Split a tab/space/newline-delimited line into fields.
static std::vector<std::string> split_line(const std::string & line)
{
    std::vector<std::string> fields;
    std::string field;
    for (char c : line) {
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            if (!field.empty()) {
                fields.push_back(field);
                field.clear();
            }
        } else {
            field += c;
        }
    }
    if (!field.empty())
        fields.push_back(field);
    return fields;
}

// ---------------------------------------------------------------------------
// Gene data loader
// ---------------------------------------------------------------------------

// Reads the gene-list file, inserts into hash table, then for each chromosome
// sorts by start position and builds the gene_order_list (start, end,
// running_max_right) triples used by the SNP-in-gene filter.
static void structure_gene_data(HashTable & ht,
                                const std::string & filename,
                                bool has_header,
                                unsigned int extend_len)
{
    std::ifstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open gene list file: " + filename);

    std::string line;

    if (has_header) {
        if (!std::getline(f, line))
            throw std::runtime_error("Gene file is empty (expected header).");
    }

    while (std::getline(f, line)) {
        if (line.empty()) continue;

        auto fields = split_line(line);
        // Expected columns: chrom, start, end, name  (0-based)
        if (fields.size() < 4) continue;

        std::string key = chrom_key(fields[0]);
        GeneNode g;
        g.gene_start = std::stoul(fields[1]);
        g.gene_end   = std::stoul(fields[2]);
        g.gene_name  = fields[3];

        ht[key].genes.push_back(std::move(g));
    }

    // Sort each chromosome's genes by start, then build the order list.
    for (auto & [key, bucket] : ht) {
        std::sort(bucket.genes.begin(), bucket.genes.end(),
                  [](const GeneNode & a, const GeneNode & b){
                      return a.gene_start < b.gene_start;
                  });

        unsigned long most_right = 0;
        bucket.gene_order_list.reserve(bucket.genes.size() * 3);

        for (const auto & g : bucket.genes) {
            unsigned long gs = (extend_len <= g.gene_start)
                               ? g.gene_start - extend_len : 1UL;
            unsigned long ge = g.gene_end + extend_len;
            if (ge > most_right) most_right = ge;
            bucket.gene_order_list.push_back(gs);
            bucket.gene_order_list.push_back(ge);
            bucket.gene_order_list.push_back(most_right);
        }
    }
}

// ---------------------------------------------------------------------------
// SNP-in-gene pre-filter  (identical logic to original judge_snp_in_gene)
// ---------------------------------------------------------------------------

static bool snp_in_gene_range(unsigned long snp_pos,
                               const std::vector<unsigned long> & order_list)
{
    const std::size_t len = order_list.size();   // already * 3
    unsigned long most_right = 0;

    for (std::size_t i = 0; i < len; i += 3) {
        if (snp_pos < order_list[i]) {
            return snp_pos > most_right ? false : true;
        }
        most_right = order_list[i + 2];
    }
    return snp_pos <= most_right;
}

// ---------------------------------------------------------------------------
// SNP data loader
// ---------------------------------------------------------------------------

static void structure_snp_data(HashTable & ht,
                                const std::string & filename,
                                bool has_header,
                                double max_af,
                                unsigned int min_sample)
{
    std::ifstream f(filename);
    if (!f)
        throw std::runtime_error("Cannot open SNP list file: " + filename);

    std::string line;

    if (has_header) {
        if (!std::getline(f, line))
            throw std::runtime_error("SNP file is empty (expected header).");
    }

    while (std::getline(f, line)) {
        if (line.empty()) continue;

        auto fields = split_line(line);
        // Column indices (0-based) matching the original:
        //  0: chrom,  2: pos,  5: sample_num,  6: af,  12: p_value
        if (fields.size() < 13) continue;

        std::string key        = chrom_key(fields[0]);
        unsigned long snp_pos  = std::stoul(fields[2]);
        unsigned long sample_n = std::stoul(fields[5]);
        double af              = std::stod(fields[6]);
        double p_value         = std::stod(fields[12]);

        af       = (af < 0.5) ? af : 1.0 - af;
        sample_n = static_cast<unsigned long>(af * sample_n * 2);

        // Apply filters; also skip chromosomes that have no gene data.
        auto it = ht.find(key);
        if (it == ht.end()) continue;
        ChromBucket & bucket = it->second;
        if (bucket.genes.empty()) continue;

        if ((af - max_af) > 1e-15) continue;
        if (sample_n < static_cast<unsigned long>(min_sample)) continue;
        if (!snp_in_gene_range(snp_pos, bucket.gene_order_list)) continue;

        SnpNode s;
        s.snp_pos = snp_pos;
        s.weight  = std::pow(BETA_DENSITY(af, 1, 25) /
                             BETA_DENSITY(af, 0.5, 0.5), 2);
        s.p_value = p_value;
        bucket.snps.push_back(s);
    }
}

// ---------------------------------------------------------------------------
// SNP → gene mapping  (replaces thread_worker_1)
// ---------------------------------------------------------------------------

// For each SNP, find every gene whose extended interval covers it and push a
// pointer into that gene's snp_ptrs list.  Genes are sorted by start; we use
// the gene_order_list triples for the same backtracking logic as the original.
static void map_snps_to_genes(ChromBucket & bucket)
{
    const auto & order = bucket.gene_order_list;
    const std::size_t list_len = order.size();   // already n_genes * 3

    for (SnpNode & snp : bucket.snps) {
        unsigned long snp_pos = snp.snp_pos;
        bool past_all = true;
        std::size_t gene_index = 0;   // tracks position in order list

        for (std::size_t i = 0; i < list_len; i += 3) {
            if (snp_pos < order[i]) {
                past_all = false;
                // Backtrack: assign to any overlapping earlier genes.
                for (std::size_t gi = gene_index; gi > 0; gi -= 3) {
                    if (snp_pos <= order[gi - 1]) {          // within end
                        if (snp_pos >= order[gi - 3]) {      // within start (gi-3 = start of that triple)
                            std::size_t gene_idx = (gi - 3) / 3;
                            bucket.genes[gene_idx].snp_ptrs.push_back(&snp);
                        }
                    } else {
                        break;
                    }
                }
                break;
            }
            gene_index += 3;
        }

        if (past_all) {
            // SNP is past all gene starts — backtrack from the end.
            for (std::size_t gi = list_len; gi > 0; gi -= 3) {
                if (snp_pos <= order[gi - 1]) {              // within end
                    if (snp_pos >= order[gi - 3]) {          // within start
                        std::size_t gene_idx = (gi - 3) / 3;
                        bucket.genes[gene_idx].snp_ptrs.push_back(&snp);
                    }
                } else {
                    break;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Result output
// ---------------------------------------------------------------------------

static void print_res_file(const HashTable & ht, const std::string & out_name)
{
    std::ofstream f(out_name);
    if (!f)
        throw std::runtime_error("Cannot open output file: " + out_name);

    f << "CHR\tGENE\tSTART\tEND\tSNP_NUM\tP_ACAT\n";

    // Iterate in a deterministic order (sort chromosome keys).
    std::vector<std::string> keys;
    keys.reserve(ht.size());
    for (const auto & [k, _] : ht) keys.push_back(k);
    std::sort(keys.begin(), keys.end());

    for (const auto & key : keys) {
        const ChromBucket & bucket = ht.at(key);
        if (bucket.genes.empty() || bucket.snps.empty()) continue;

        for (const auto & gene : bucket.genes) {
            if (gene.snp_ptrs.empty()) continue;

            double weight_sum  = 0.0;
            double c_value_sum = 0.0;
            unsigned long snp_amount = gene.snp_ptrs.size();

            for (const SnpNode * s : gene.snp_ptrs)
                weight_sum += s->weight;

            for (const SnpNode * s : gene.snp_ptrs) {
                double p      = s->p_value;
                double w_ave  = s->weight / weight_sum;
                if (p >= 1e-15) {
                    double tmp = (0.5 - p) * M_PI;
                    c_value_sum += w_ave * (std::sin(tmp) / std::cos(tmp));
                } else {
                    c_value_sum += (w_ave / p) * M_PI;
                }
            }

            double cauchy;
            if (std::fabs(c_value_sum) > 1.0) {
                cauchy = PCAUCHY2(c_value_sum);
                cauchy = (c_value_sum > 0) ? cauchy : 1.0 + cauchy;
            } else {
                cauchy = 1.0L - PCAUCHY(c_value_sum);
            }

            f << key << '\t'
              << gene.gene_name << '\t'
              << gene.gene_start << '\t'
              << gene.gene_end   << '\t'
              << snp_amount      << '\t'
              << std::scientific << cauchy << '\n';
        }
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

int acat_func(const std::string & gene_list_file,
              const std::string & snp_list_file,
              double              max_af,
              unsigned int        min_sample,
              unsigned int        extend_len,
              const std::string & out_f)
{
    std::cout << "args: --acat --gene-list " << gene_list_file
              << " --snp-list "  << snp_list_file
              << " --max-maf "   << max_af
              << " --min-mac "   << min_sample
              << " --wind "      << extend_len
              << " --out "       << out_f << '\n';

    try {
        HashTable ht;
        time_t t0, t1;

        std::cout << ">structure gene_list data start...\n";
        time(&t0);
        structure_gene_data(ht, gene_list_file, false, extend_len);
        time(&t1);
        std::cout << "<done. used " << difftime(t1, t0) << " second(s)\n\n";

        std::cout << ">structure snp_list data start...\n";
        time(&t0);
        structure_snp_data(ht, snp_list_file, true, max_af, min_sample);
        time(&t1);
        std::cout << "<done. used " << difftime(t1, t0) << " second(s)\n\n";

        std::cout << ">mapping snp to gene start...\n";
        time(&t0);
        for (auto & [key, bucket] : ht) {
            if (!bucket.genes.empty() && !bucket.snps.empty())
                map_snps_to_genes(bucket);
        }
        time(&t1);
        std::cout << "<done. used " << difftime(t1, t0) << " second(s)\n\n";

        std::cout << ">calculate cauchy and print results start...\n";
        time(&t0);
        print_res_file(ht, out_f);
        time(&t1);
        std::cout << "<done. used " << difftime(t1, t0) << " second(s)\n\n";

        // Memory released automatically when ht goes out of scope.
        return 0;

    } catch (const std::exception & e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}

// ---------------------------------------------------------------------------
// Optional standalone entry point
// ---------------------------------------------------------------------------

#ifdef ACAT_UTILITY
int main(int argc, char * argv[])
{
    if (argc != 3) {
        std::cerr << "error, argument number error.\n";
        return EXIT_FAILURE;
    }
    std::string  gene_file  = argv[1];
    std::string  snp_file   = argv[2];
    double       max_af     = 0.01;
    unsigned int min_sample = 10;
    unsigned int expand_len = 0;
    std::string  res_file   = "test_result.csv";

    return acat_func(gene_file, snp_file, max_af, min_sample, expand_len, res_file);
}
#endif