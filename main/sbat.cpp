/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementation of gene-based association test (GBAT) in GCTA
 *
 * 2013 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "gcta.h"
#include <set>

void gcta::sbat_read_snpAssoc(std::string snpAssoc_file, std::vector<std::string> &snp_name, std::vector<int> &snp_chr, std::vector<int> &snp_bp, std::vector<double> &snp_pval)
{
    std::ifstream in_snpAssoc(snpAssoc_file.c_str());
    if (!in_snpAssoc) LOGGER.e(0, "cannot open the file [" + snpAssoc_file + "] to read.");
    LOGGER << "\nReading SNP association results from [" + snpAssoc_file + "]." << std::endl;
    std::string str_buf;
    std::vector<std::string> vs_buf;
    std::map<std::string, int>::iterator iter;
    std::map<std::string, int> assoc_snp_map;
    int line = 0;
    while (std::getline(in_snpAssoc, str_buf)) {
        if (StrFunc::split_string(str_buf, vs_buf, " \t") != 2) LOGGER.e(0, "in line \"" + str_buf + "\".");
        iter = _snp_name_map.find(vs_buf[0]);
        if (iter == _snp_name_map.end()) continue;
        if(assoc_snp_map.find(vs_buf[0]) != assoc_snp_map.end()) continue;
        else assoc_snp_map.insert(std::pair<std::string, int>(vs_buf[0], line));
        snp_name.push_back(vs_buf[0]);
        snp_pval.push_back(atof(vs_buf[1].c_str()));
        line++;
    }
    in_snpAssoc.close();
    snp_name.erase(unique(snp_name.begin(), snp_name.end()), snp_name.end());
    LOGGER << "Association p-values of " << snp_name.size() << " SNPs have been included." << std::endl;

    update_id_map_kp(snp_name, _snp_name_map, _include);
    std::vector<std::string> snp_name_buf(snp_name);
    std::vector<double> snp_pval_buf(snp_pval);
    snp_name.clear();
    snp_pval.clear();
    snp_name.resize(_include.size());
    snp_pval.resize(_include.size());
    int i = 0;
    std::map<std::string, int> snp_name_buf_map;
    for (i = 0; i < snp_name_buf.size(); i++) snp_name_buf_map.insert(std::pair<std::string,int>(snp_name_buf[i], i));
    #pragma omp parallel for
    for (i = 0; i < _include.size(); i++) {
        std::map<std::string, int>::iterator iter = snp_name_buf_map.find(_snp_name[_include[i]]);
        snp_name[i] = snp_name_buf[iter->second];
        snp_pval[i] = snp_pval_buf[iter->second];
    }
    snp_chr.resize(_include.size());
    snp_bp.resize(_include.size());
    #pragma omp parallel for
    for (i = 0; i < _include.size(); i++) {
        snp_chr[i] = _chr[_include[i]];
        snp_bp[i] = _bp[_include[i]];
    }
    if (_include.size() < 1) LOGGER.e(0, "no SNP is included in the analysis.");
    else if (_chr[_include[0]] < 1) LOGGER.e(0, "chromosome information is missing.");
    else if (_bp[_include[0]] < 1) LOGGER.e(0, "bp information is missing.");
}

void gcta::sbat_read_geneAnno(std::string gAnno_file, std::vector<std::string> &gene_name, std::vector<int> &gene_chr, std::vector<int> &gene_bp1, std::vector<int> &gene_bp2) {
    std::ifstream in_gAnno(gAnno_file.c_str());
    if (!in_gAnno) LOGGER.e(0, "cannot open the file [" + gAnno_file + "] to read.");
    LOGGER << "Reading physical positions of the genes from [" + gAnno_file + "]." << std::endl;
    std::string str_buf;
    std::vector<std::string> vs_buf;
    while (std::getline(in_gAnno, str_buf)) {
        if (StrFunc::split_string(str_buf, vs_buf) != 4) LOGGER.e(0, "in line \"" + str_buf + "\".");
        gene_chr.push_back(atoi(vs_buf[0].c_str()));
        gene_bp1.push_back(atoi(vs_buf[1].c_str()));
        gene_bp2.push_back(atoi(vs_buf[2].c_str()));
        gene_name.push_back(vs_buf[3]);
    }
    in_gAnno.close();
    LOGGER << "Physical positions of " << gene_name.size() << " genes have been include." << std::endl;
}

void gcta::sbat_gene(std::string sAssoc_file, std::string gAnno_file, int wind, double sbat_ld_cutoff, bool sbat_write_snpset, bool GC, double GC_val)
{
    int i = 0, j = 0;
    int snp_count;

    // read SNP association results
    // std::vector<std::string> snp_name;
    // std::vector<int> snp_chr, snp_bp;
    // std::vector<double> snp_pval;
    // sbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
    // std::vector<double> snp_chisq(snp_pval.size());
    // for (i = 0; i < snp_pval.size(); i++) snp_chisq[i] = StatFunc::qchisq(snp_pval[i], 1);

    // read SNP association results
    std::vector<std::string> snp_name;
    std::vector<int> snp_chr, snp_bp;
    std::vector<double> snp_pval;
    init_massoc(sAssoc_file, GC, GC_val);
    int snp_num = _include.size();
    // re-calculate chi-square 
    std::vector<double> snp_chisq(snp_num);
    snp_pval.resize(snp_num);
    snp_name.resize(snp_num);
    snp_chr.resize(snp_num);
    snp_bp.resize(snp_num);
    for (i = 0; i < snp_num; i++) {
        snp_name[i] = _snp_name[_include[i]];
        snp_chr[i]  = _chr[_include[i]];
        snp_bp[i] =  _bp[_include[i]];
        snp_pval[i] = _pval[i];
        snp_chisq[i] = StatFunc::qchisq(_pval[i], 1);
    }
    
    // get start and end of chr
    // int snp_num = snp_name.size();
    std::map<int, std::string> chr_begin_snp, chr_end_snp;
    chr_begin_snp.insert(std::pair<int, std::string>(snp_chr[0], snp_name[0]));
    for (i = 1; i < snp_num; i++) {
        if (snp_chr[i] != snp_chr[i - 1]) {
            chr_begin_snp.insert(std::pair<int, std::string>(snp_chr[i], snp_name[i]));
            chr_end_snp.insert(std::pair<int, std::string>(snp_chr[i - 1], snp_name[i - 1]));
        }
    }
    chr_end_snp.insert(std::pair<int, std::string>(snp_chr[snp_num - 1], snp_name[snp_num - 1]));
    
    // read gene list
    std::vector<std::string> gene_name;
    std::vector<int> gene_chr, gene_bp1, gene_bp2;
    sbat_read_geneAnno(gAnno_file, gene_name, gene_chr, gene_bp1, gene_bp2);

    // std::map genes to SNPs
    LOGGER << "Mapping the physical positions of genes to SNP data (gene boundaries: " << wind / 1000 << "Kb away from UTRs) ..." << std::endl;
    int gene_num = gene_name.size();
    std::vector<std::string> gene2snp_1(gene_num), gene2snp_2(gene_num);
    std::vector<locus_bp>::iterator iter;
    std::map<int, std::string>::iterator chr_iter;
    std::vector<locus_bp> snp_vec;
    for (i = 0; i < snp_num; i++) snp_vec.push_back(locus_bp(snp_name[i], snp_chr[i], snp_bp[i]));
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp1[i] - wind));
        if (iter != snp_vec.end()) gene2snp_1[i] = iter->locus_name;
        else gene2snp_1[i] = "NA";
    }
    #pragma omp parallel for private(iter, chr_iter)
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] == "NA") {
            gene2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snp_vec.begin(), snp_vec.end(), locus_bp(gene_name[i], gene_chr[i], gene_bp2[i] + wind));
        if (iter != snp_vec.end()){
            if (iter->bp ==  gene_bp2[i] + wind) gene2snp_2[i] = iter->locus_name;
            else {
                if(iter!=snp_vec.begin()){
                    iter--;
                    gene2snp_2[i] = iter->locus_name;
                }
                else gene2snp_2[i] = "NA";
            }
        }
        else {
            chr_iter = chr_end_snp.find(gene_chr[i]);
            if (chr_iter == chr_end_snp.end()) gene2snp_2[i] = "NA";
            else gene2snp_2[i] = chr_iter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < gene_num; i++) {
        if (gene2snp_1[i] != "NA" && gene2snp_2[i] != "NA") mapped++;
    }
    if (mapped < 1) LOGGER.e(0, "no gene can be mapped to the SNP data. Please check the input data regarding chromosome and bp.");
    else LOGGER << mapped << " genes have been mapped to SNP data." << std::endl;

    // run gene-based test
    if (_mu.empty()) calcu_mu();
    LOGGER << "\nRunning fastBAT analysis for genes ..." << std::endl;
    if (sbat_ld_cutoff < 1) LOGGER << "Pruning SNPs with LD rsq cutoff = " << sbat_ld_cutoff*sbat_ld_cutoff  << std::endl;
    std::vector<double> gene_pval(gene_num), chisq_o(gene_num), min_snp_pval(gene_num),eigenval_fastbat(gene_num);
    std::vector<std::string> min_snp_name(gene_num);
    std::vector<int> snp_num_in_gene(gene_num);
    std::map<std::string, int>::iterator iter1, iter2;
    std::map<std::string, int> snp_name_map;
    std::string rgoodsnpfile = _out + ".gene.snpset";
    std::ofstream rogoodsnp;
    
    if (sbat_write_snpset) rogoodsnp.open(rgoodsnpfile.c_str());
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(std::pair<std::string,int>(snp_name[i], i));
    for (i = 0; i < gene_num; i++) {
        iter1 = snp_name_map.find(gene2snp_1[i]);
        iter2 = snp_name_map.find(gene2snp_2[i]);
        bool skip = false;
        if (iter1 == snp_name_map.end() || iter2 == snp_name_map.end() || iter1->second >= iter2->second) skip = true;
        snp_num_in_gene[i] = iter2->second - iter1->second + 1;
        if(!skip && snp_num_in_gene[i] > 10000){
            LOGGER<<"Warning: Too many SNPs in the gene region ["<<gene_name[i]<<"]. Maximum limit is 10000. This gene is ignored in the analysis."<<std::endl;
            skip = true;  
        } 
        if(skip){
            gene_pval[i] = 2.0;
            snp_num_in_gene[i] = 0;
            continue;
        }
        chisq_o[i] = 0;
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = iter1->second; j <= iter2->second; j++) {
            if (min_snp_pval[i] > snp_pval[j]) { 
                min_snp_pval[i] = snp_pval[j];
                min_snp_name[i] = snp_name[j]; //keep minimum value - regardless of whether SNP removed by LD pruning
            }
            chisq_o[i] += snp_chisq[j];
        }
        if(snp_num_in_gene[i] == 1) gene_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            std::vector<int> snp_indx;
            for (j = iter1->second; j <= iter2->second; j++) snp_indx.push_back(j);            
            snp_count=snp_num_in_gene[i];
            Eigen::VectorXd eigenval;
            std::vector<int> sub_indx;
            sbat_calcu_lambda(snp_indx, eigenval, snp_count, sbat_ld_cutoff, sub_indx);
            //recalculate chisq value from low correlation snp subset
            // eigenval_fastbat[i] = eigenval;
            if (sbat_ld_cutoff < 1) {
                chisq_o[i] = 0;
                for (j = 0; j < sub_indx.size(); j++) chisq_o[i] += snp_chisq[snp_indx[sub_indx[j]]];
            } 
            snp_num_in_gene[i] = snp_count;
            if (snp_count==1 && chisq_o[i] ==0) gene_pval[i] = 1;
            else gene_pval[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);

            if (sbat_write_snpset) {
                rogoodsnp << gene_name[i] << std::endl;
                for (int k = 0; k < sub_indx.size(); k++) rogoodsnp << snp_name[(iter1->second)+sub_indx[k]] << std::endl;
                rogoodsnp << "END" << std::endl << std::endl;
            }

        }

        if((i + 1) % 100 == 0 || (i + 1) == gene_num) LOGGER << i + 1 << " of " << gene_num << " genes.\r";
    }

    std::string filename = _out + ".gene.fastbat";
    LOGGER << "\nSaving the results of the fastBAT analysis to [" + filename + "] ..." << std::endl;
    std::ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Gene\tChr\tStart\tEnd\tNo.SNPs\tSNP_start\tSNP_end\tChisq(Obs)\tPvalue\tTopSNP.Pvalue\tTopSNP" << std::endl;
    for (i = 0; i < gene_num; i++) {
        if(gene_pval[i]>1.5) continue;
        ofile << gene_name[i] << "\t" << gene_chr[i] << "\t" << gene_bp1[i] << "\t" << gene_bp2[i] << "\t";
        ofile << snp_num_in_gene[i] << "\t" << gene2snp_1[i] << "\t" << gene2snp_2[i] << "\t" << chisq_o[i]; // << "\t" << eigenval_fastbat[i];
        ofile << "\t" << gene_pval[i] << "\t" << min_snp_pval[i] << "\t" << min_snp_name[i] << std::endl;
        //else ofile << "0\tNA\tNA\tNA\tNA" << std::endl;
    }
    ofile.close();
    if (sbat_write_snpset) {
        LOGGER << "The SNP sets have been saved in file [" << rgoodsnpfile << "]." << std::endl;
        rogoodsnp.close();
    }
}


void gcta::sbat_read_snpset(std::string snpset_file, std::vector<std::string> &set_name, std::vector< std::vector<std::string> > &snpset)
{
    std::ifstream in_snpset(snpset_file.c_str());
    if (!in_snpset) LOGGER.e(0, "cannot open the file [" + snpset_file + "] to read.");
    LOGGER << "\nReading SNP sets from [" + snpset_file + "]." << std::endl;
    std::string str_buf;
    std::vector<std::string> vs_buf, snpset_buf, snp_name;
    // add std::set std::vector to remove duplicate items;
    std::set<std::string> snpUniqSet_buf;
    int i = 0;
    while (in_snpset>>str_buf) {
        if(str_buf!="END" && str_buf!="end") vs_buf.push_back(str_buf);
        else{
            if(vs_buf.empty()) continue;
            set_name.push_back(vs_buf[0]);
            snpUniqSet_buf.clear();
            for(i = 1; i < vs_buf.size(); i++){
                if (_snp_name_map.find(vs_buf[i]) != _snp_name_map.end()){
                    // if there is duplicated snp, ignore it
                    if(snpUniqSet_buf.find(vs_buf[i]) != snpUniqSet_buf.end()) continue;
                    snpUniqSet_buf.insert(vs_buf[i]);

                    snpset_buf.push_back(vs_buf[i]);
                    snp_name.push_back(vs_buf[i]);
                    
                }
            }
            vs_buf.clear();
            snpset.push_back(snpset_buf);
            snpset_buf.clear();
        }
    }
    in_snpset.close();
    snp_name.erase(unique(snp_name.begin(), snp_name.end()), snp_name.end());
    update_id_map_kp(snp_name, _snp_name_map, _include);
    // cout  << "snpset[0] size: " << snpset[0].size() << std::endl;
    LOGGER << snp_name.size() << " SNPs in " << snpset.size() << " sets have been included." << std::endl;
}


void gcta::sbat(std::string sAssoc_file, std::string snpset_file, double sbat_ld_cutoff, bool sbat_write_snpset,bool GC, double GC_val)

{
    int i = 0, j = 0;
    int snp_count;

    // read SNP std::set file
    std::vector<std::string> set_name;
    std::vector< std::vector<std::string> > snpset;
    sbat_read_snpset(snpset_file, set_name, snpset);
    int set_num = set_name.size();

    // read SNP association results
    std::vector<std::string> snp_name;
    std::vector<int> snp_chr, snp_bp;
    std::vector<double> snp_pval;
    init_massoc(sAssoc_file, GC, GC_val);
    int snp_num = _include.size();
    // re-calculate chi-square 
    std::vector<double> snp_chisq(snp_num);
    snp_pval.resize(snp_num);
    snp_name.resize(snp_num);
    snp_chr.resize(snp_num);
    snp_bp.resize(snp_num);
    for (i = 0; i < snp_num; i++) {
        snp_name[i] = _snp_name[_include[i]];
        snp_chr[i]  = _chr[_include[i]];
        snp_bp[i] =  _bp[_include[i]];
        snp_pval[i] = _pval[i];
        snp_chisq[i] = StatFunc::qchisq(_pval[i], 1);
    }
    // run gene-based test
    if (_mu.empty()) calcu_mu();
    LOGGER << "\nRunning fastBAT analysis ..." << std::endl;
    if (sbat_ld_cutoff < 1) LOGGER << "Pruning SNPs with maximum LD cutoff " << sbat_ld_cutoff  << std::endl;
    std::vector<double> set_pval(set_num), chisq_o(set_num), min_snp_pval(set_num);
    std::vector<std::string> min_snp_name(set_num);
    std::vector<int> snp_num_in_set(set_num);
    std::map<std::string, int>::iterator iter;
    std::map<std::string, int> snp_name_map;

    std::string rgoodsnpfile = _out + ".snpset";
    std::ofstream rogoodsnp;
    if (sbat_write_snpset) rogoodsnp.open(rgoodsnpfile.c_str());
 
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(std::pair<std::string,int>(snp_name[i], i));
    for (i = 0; i < set_num; i++) {
        bool skip = false;
        if(snpset[i].size() < 1) skip = true;
        std::vector<int> snp_indx;
        for(j = 0; j < snpset[i].size(); j++){
            iter = snp_name_map.find(snpset[i][j]);
            if(iter!=snp_name_map.end()) snp_indx.push_back(iter->second);
        }
        snp_num_in_set[i] = snp_indx.size();
        if(!skip && snp_num_in_set[i] > 20000){
            LOGGER<<"Warning: Too many SNPs in the std::set ["<<set_name[i]<<"]. Maximum limit is 20000. This gene is ignored in the analysis."<<std::endl;
            skip = true;  
        } 
        if(skip){
            set_pval[i] = 2.0;
            snp_num_in_set[i] = 0;
            continue;
        }
        chisq_o[i] = 0;
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = 0; j < snp_indx.size(); j++) {
            if (min_snp_pval[i] > snp_pval[snp_indx[j]]) {
                min_snp_pval[i] = snp_pval[snp_indx[j]];
                min_snp_name[i] = snp_name[snp_indx[j]];
            }
            chisq_o[i] += snp_chisq[snp_indx[j]];
        }
        if(snp_num_in_set[i] == 1) set_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            snp_count=snp_num_in_set[i];
            Eigen::VectorXd eigenval;
            std::vector<int> sub_indx;
            sbat_calcu_lambda(snp_indx, eigenval, snp_count, sbat_ld_cutoff, sub_indx);

            //recalculate chisq value from low correlation snp subset
            if (sbat_ld_cutoff < 1) {
                chisq_o[i] = 0;
                for (j = 0; j < sub_indx.size(); j++) chisq_o[i] += snp_chisq[snp_indx[sub_indx[j]]];
            }
            snp_num_in_set[i] = snp_count;
            if (snp_count==1 && chisq_o[i] ==0) set_pval[i] = 1;
            else set_pval[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);

            if (sbat_write_snpset) {
                rogoodsnp << set_name[i] << std::endl;
                for (int k = 0; k < sub_indx.size(); k++) rogoodsnp << snp_name[snp_indx[sub_indx[k]]] << std::endl;
                rogoodsnp << "END" << std::endl << std::endl;
            }

        }

        if((i + 1) % 100 == 0 || (i + 1) == set_num) LOGGER << i + 1 << " of " << set_num << " sets.\r";
    }

    std::string filename = _out + ".fastbat";
    LOGGER << "\nSaving the results of the fastBAT analysis to [" + filename + "] ..." << std::endl;
    std::ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Set\tNo.SNPs\tChisq(Obs)\tPvalue\tTopSNP.Pvalue\tTopSNP" << std::endl;
    for (i = 0; i < set_num; i++) {
        if(set_pval[i]>1.5) continue;
        ofile << set_name[i] << "\t" << snp_num_in_set[i] << "\t" << chisq_o[i] << "\t";
        ofile << set_pval[i] << "\t" << min_snp_pval[i] << "\t" << min_snp_name[i] << std::endl;
    }
    ofile.close();
    if (sbat_write_snpset) {
        LOGGER << "The SNP sets have been saved in file [" << rgoodsnpfile << "]." << std::endl;
        rogoodsnp.close();
    }
}


void gcta::sbat_seg(std::string sAssoc_file, int seg_size, double sbat_ld_cutoff, bool sbat_write_snpset,bool GC, double GC_val)
{
    int i = 0, j = 0;
    int snp_count;

    // read SNP association results

    std::vector<std::string> snp_name;
    std::vector<int> snp_chr, snp_bp;
    std::vector<double> snp_pval;
    init_massoc(sAssoc_file, GC, GC_val);
    int snp_num = _include.size();
    // re-calculate chi-square 
    std::vector<double> snp_chisq(snp_num);
    snp_pval.resize(snp_num);
    snp_name.resize(snp_num);
    snp_chr.resize(snp_num);
    snp_bp.resize(snp_num);
    for (i = 0; i < snp_num; i++) {
        snp_name[i] = _snp_name[_include[i]];
        snp_chr[i]  = _chr[_include[i]];
        snp_bp[i] =  _bp[_include[i]];
        snp_pval[i] = _pval[i];
        snp_chisq[i] = StatFunc::qchisq(_pval[i], 1);
    }

    // run gene-based test
    if (_mu.empty()) calcu_mu();
    LOGGER << "\nRunning fastBAT analysis at genomic segments with a length of " << seg_size/1000 << "Kb ..." << std::endl;
    if (sbat_ld_cutoff < 1) LOGGER << "Pruning SNPs with maximum LD cutoff " << sbat_ld_cutoff  << std::endl;
    std::vector< std::vector<int> > snp_set_indx;
    std::vector<int> set_chr, set_start_bp, set_end_bp;
    get_sbat_seg_blk(seg_size, snp_set_indx, set_chr, set_start_bp, set_end_bp);
    int set_num = snp_set_indx.size();
    std::vector<double> set_pval(set_num), chisq_o(set_num), min_snp_pval(set_num);
    std::vector<std::string> min_snp_name(set_num);
    std::vector<int> snp_num_in_set(set_num);

    std::string rgoodsnpfile = _out + ".seg.snpset";
    std::ofstream rogoodsnp;
    if (sbat_write_snpset) rogoodsnp.open(rgoodsnpfile.c_str());
 
    for (i = 0; i < set_num; i++) {
        bool skip = false;
        std::vector<int> snp_indx = snp_set_indx[i];
        if(snp_indx.size() < 1) skip = true;
        snp_num_in_set[i] = snp_indx.size();
        if(!skip && snp_num_in_set[i] > 20000){
            LOGGER<<"Warning: Too many SNPs in the std::set on [chr" << set_chr[i] << ":" << set_start_bp[i] << "-" << set_end_bp[i] << "]. Maximum limit is 20000. This gene is ignored in the analysis."<<std::endl;
            skip = true;  
        } 
        if(skip){
            set_pval[i] = 2.0;
            snp_num_in_set[i] = 0;
            continue;
        }
        chisq_o[i] = 0; 
        min_snp_pval[i]=2;
        min_snp_name[i]="na";
        for (j = 0; j < snp_indx.size(); j++) {
            if (min_snp_pval[i] > snp_pval[snp_indx[j]]) {
                min_snp_pval[i] = snp_pval[snp_indx[j]];
                min_snp_name[i] = snp_name[snp_indx[j]];
            }
            chisq_o[i] += snp_chisq[snp_indx[j]];
        }
        if(snp_num_in_set[i] == 1) set_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            snp_count=snp_num_in_set[i];
            Eigen::VectorXd eigenval;
            std::vector<int> sub_indx;
            sbat_calcu_lambda(snp_indx, eigenval, snp_count, sbat_ld_cutoff, sub_indx);
            //recalculate chisq value from low correlation snp subset
            if (sbat_ld_cutoff < 1) {
                chisq_o[i] = 0;
                for (j = 0; j < sub_indx.size(); j++) chisq_o[i] += snp_chisq[snp_indx[sub_indx[j]]]; 
            }
            snp_num_in_set[i] = snp_count;
            if (snp_count==1 && chisq_o[i] ==0) set_pval[i] = 1;
            else set_pval[i] = StatFunc::pchisqsum(chisq_o[i], eigenval);

           if (sbat_write_snpset) {
                rogoodsnp << "SEG" << i << ":" << set_start_bp[i] << "-" << set_end_bp[i] << std::endl;
                for (int k = 0; k < sub_indx.size(); k++) rogoodsnp << snp_name[snp_indx[sub_indx[k]]] << std::endl;
                rogoodsnp << "END" << std::endl << std::endl;
            }

        }

        if((i + 1) % 100 == 0 || (i + 1) == set_num) LOGGER << i + 1 << " of " << set_num << " sets.\r";
    }

    std::string filename = _out + ".seg.fastbat";
    LOGGER << "\nSaving the results of the segment-based fastBAT analysis to [" + filename + "] ..." << std::endl;
    std::ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Chr\tStart\tEnd\tNo.SNPs\tChisq(Obs)\tPvalue\tTopSNP.Pvalue\tTopSNP" << std::endl;
    for (i = 0; i < set_num; i++) {
        if(set_pval[i]>1.5) continue;
        ofile << set_chr[i] << "\t" << set_start_bp[i] << "\t"<< set_end_bp[i] << "\t";
        ofile << snp_num_in_set[i] << "\t" << chisq_o[i] << "\t" << set_pval[i] << "\t";
        ofile << min_snp_pval[i] << "\t" << min_snp_name[i] << std::endl;
 
    }
    ofile.close();
    if (sbat_write_snpset) {
        LOGGER << "The SNP sets have been saved in file [" << rgoodsnpfile << "]." << std::endl;
        rogoodsnp.close();
    }
}


void gcta::get_sbat_seg_blk(int seg_size, std::vector< std::vector<int> > &snp_set_indx, std::vector<int> &set_chr, std::vector<int> &set_start_bp, std::vector<int> &set_end_bp)
{
    int i = 0, j = 0, k = 0, m = _include.size();

    std::vector<int> brk_pnt;
    brk_pnt.push_back(0);
    for (i = 1, j = 0; i < m; i++) {
        if (i == (m - 1)) brk_pnt.push_back(m - 1);
        else if (_chr[_include[i]] != _chr[_include[brk_pnt[j]]]) {
            brk_pnt.push_back(i - 1);
            j++;
            brk_pnt.push_back(i);
            j++;
        }
        else if (_bp[_include[i]] - _bp[_include[brk_pnt[j]]] > seg_size) {
            brk_pnt.push_back(i - 1);
            j++;
            brk_pnt.push_back(i);
            j++;
        }
    }

    snp_set_indx.clear();
    set_chr.clear();
    set_start_bp.clear();
    set_end_bp.clear();
    for (i = 0; i < brk_pnt.size() - 1; i++) {
        int size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if(size < 3 && (i%2 != 0)) continue;
        std::vector<int> snp_indx(size);
        for (j = brk_pnt[i], k = 0; j <= brk_pnt[i + 1]; j++, k++){
            snp_indx[k] = j;  
        } 
        snp_set_indx.push_back(snp_indx);
        set_chr.push_back(_chr[_include[brk_pnt[i]]]);
        set_start_bp.push_back(_bp[_include[brk_pnt[i]]]);
        set_end_bp.push_back(_bp[_include[brk_pnt[i + 1]]]);
    }
}

void gcta::sbat_calcu_lambda(std::vector<int> &snp_indx, Eigen::VectorXd &eigenval, int &snp_count, double sbat_ld_cutoff, std::vector<int> &sub_indx)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = snp_indx.size();

    Eigen::MatrixXf X;
    make_XMat_subset(X, snp_indx, false);
    std::vector<int> rm_ID1;
    double R_cutoff = sbat_ld_cutoff;
    int qi = 0; //alternate index

    Eigen::VectorXd sumsq_x(m);
    for (j = 0; j < m; j++) sumsq_x[j] = X.col(j).dot(X.col(j));

    Eigen::MatrixXf C = X.transpose() * X;
    X.resize(0,0);
    #pragma omp parallel for private(j)
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            double d_buf = sqrt(sumsq_x[i] * sumsq_x[j]);
            if(d_buf>0.0) C(i,j) /= d_buf;
            else C(i,j) = 0.0;
        }
    }

    if (sbat_ld_cutoff < 1) rm_cor_sbat(C, R_cutoff, m, rm_ID1);
        //Create new index
        for (int i=0 ; i<m ; i++) {
            if (rm_ID1.size() == 0) sub_indx.push_back(i);
            else {
                if (rm_ID1[qi] == i) qi++; //Skip removed snp
                else sub_indx.push_back(i);
            }
        }
        snp_count = sub_indx.size();
        if (sub_indx.size() < C.size()) { //Build new matrix
            Eigen::MatrixXf D(sub_indx.size(),sub_indx.size());
            for (i = 0 ; i < sub_indx.size() ; i++) {
               for (j = 0 ; j < sub_indx.size() ; j++) {
                   D(i,j) = C(sub_indx[i],sub_indx[j]);
               }
            }
            C = D; 
        }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> saes(C);
    eigenval = saes.eigenvalues().cast<double>();
}

void gcta::rm_cor_sbat(Eigen::MatrixXf &R, double R_cutoff, int m, std::vector<int> &rm_ID1) {
    //Modified version of rm_cor_indi from grm.cpp
    
    int i = 0, j = 0, i_buf = 0;
    std::vector<int> rm_ID2;

    //float tmpr = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < i; j++) {
            if (fabs(R(i,j)) > R_cutoff ) { 
                rm_ID1.push_back(i);
                rm_ID2.push_back(j);
            }
        }
    }

    // count the number of appearance of each "position" in the std::vector, which involves a few steps
    std::vector<int> rm_uni_ID(rm_ID1);
    rm_uni_ID.insert(rm_uni_ID.end(), rm_ID2.begin(), rm_ID2.end());
    std::stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    std::map<int, int> rm_uni_ID_count;
    for (i = 0; i < rm_uni_ID.size(); i++) {
        i_buf = count(rm_ID1.begin(), rm_ID1.end(), rm_uni_ID[i]) + count(rm_ID2.begin(), rm_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(std::pair<int, int>(rm_uni_ID[i], i_buf));
    }

    // swapping
    std::map<int, int>::iterator iter1, iter2;
    for (i = 0; i < rm_ID1.size(); i++) {
        iter1 = rm_uni_ID_count.find(rm_ID1[i]);
        iter2 = rm_uni_ID_count.find(rm_ID2[i]);
        if (iter1->second < iter2->second) {
            i_buf = rm_ID1[i];
            rm_ID1[i] = rm_ID2[i];
            rm_ID2[i] = i_buf;
        }
    }
    std::stable_sort(rm_ID1.begin(), rm_ID1.end());
    rm_ID1.erase(unique(rm_ID1.begin(), rm_ID1.end()), rm_ID1.end());
}
