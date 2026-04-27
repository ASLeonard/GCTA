#ifndef ACAT_HEAD
#define ACAT_HEAD

#include <string>

int acat_func(const std::string & gene_list_file,
              const std::string & snp_list_file,
              double              max_af,
              unsigned int        min_sample,
              unsigned int        extend_len,
              const std::string & out_f);

#endif