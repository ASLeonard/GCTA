/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementation of gene-based association test (GBAT) in GCTA
 *
 * 2013 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::gbat_read_snpAssoc(std::string snpAssoc_file, std::vector<std::string> &snp_name, std::vector<int> &snp_chr, std::vector<int> &snp_bp, std::vector<double> &snp_pval)
{
    std::ifstream in_snpAssoc(snpAssoc_file.c_str());
    if (!in_snpAssoc) LOGGER.e(0, "cannot open the file [" + snpAssoc_file + "] to read.");
    LOGGER << "\nReading SNP association results from [" + snpAssoc_file + "]." << std::endl;
    std::string str_buf;
    std::vector<std::string> vs_buf;
    LOGGER << "Reading association p-values from [" << snpAssoc_file << "]." << std::endl;
    while (std::getline(in_snpAssoc, str_buf)) {
        if (StrFunc::split_string(str_buf, vs_buf) != 2) LOGGER.e(0, "in line \"" + str_buf + "\".");
        auto iter = _snp_name_map.find(vs_buf[0]);
        if (iter == _snp_name_map.end()) continue;
        snp_name.push_back(vs_buf[0]);
        snp_pval.push_back(atof(vs_buf[1].c_str()));
    }
    in_snpAssoc.close();
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

void gcta::gbat_read_geneAnno(std::string gAnno_file, std::vector<std::string> &gene_name, std::vector<int> &gene_chr, std::vector<int> &gene_bp1, std::vector<int> &gene_bp2) {
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
    LOGGER << "Physical positions of " << gene_name.size() << " genes have been included." << std::endl;
}

void gcta::gbat_calcu_ld(Eigen::MatrixXf &X, eigenVector &sumsq_x, int snp1_indx, int snp2_indx, Eigen::MatrixXf &C)
{
    int i = 0, j = 0;
    int size = snp2_indx - snp1_indx + 1;

    C.resize(0, 0);
    if (size == 1) {
        C.resize(1, 1);
        C(1, 1) = 1.0;
        return;
    }

    //MatrixXf X_sub = X.block(0,snp1_indx,_keep.size(),size);
    C = X.block(0,snp1_indx,_keep.size(),size).transpose() * X.block(0,snp1_indx,_keep.size(),size);
    #pragma omp parallel for private(j)
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            double d_buf = sqrt(sumsq_x[snp1_indx + i] * sumsq_x[snp1_indx + j]);
            if(d_buf>0.0) C(i, j) /= d_buf;
            else C(i, j) = 0.0;
        }
    }
}

void gcta::gbat(std::string sAssoc_file, std::string gAnno_file, int wind, int simu_num)
{
    int i = 0, j = 0;

    // read SNP association results
    std::vector<std::string> snp_name;
    std::vector<int> snp_chr, snp_bp;
    std::vector<double> snp_pval;
    gbat_read_snpAssoc(sAssoc_file, snp_name, snp_chr, snp_bp, snp_pval);
    std::vector<double> snp_chisq(snp_pval.size());
    for (i = 0; i < snp_pval.size(); i++) snp_chisq[i] = StatFunc::qchisq(snp_pval[i], 1);

    // get start and end of chr
    int snp_num = snp_name.size();
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
    gbat_read_geneAnno(gAnno_file, gene_name, gene_chr, gene_bp1, gene_bp2);

    // std::map genes to SNPs
    LOGGER << "Mapping the physical positions of genes to SNP data (gene bounaries: " << wind / 1000 << "Kb away from UTRs) ..." << std::endl;

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

    // recoding genotype
    Eigen::MatrixXf X;
    eigenVector sumsq_x(_include.size());
    make_XMat(X);
    #pragma omp parallel for private(j)
    for(i = 0; i < _keep.size(); i++){
        for(j = 0; j < _include.size(); j++){
            if(X(i,j) < 1e5) X(i,j) -= _mu[_include[j]];
            else X(i,j) = 0.0;
        }
    }
    for (i = 0; i < _include.size(); i++) sumsq_x[i] = X.col(i).dot(X.col(i));

    // run gene-based test
    LOGGER << "\nRunning gene-based association test (GBAT)..." << std::endl;
    std::vector<double> gene_pval(gene_num), chisq_o(gene_num);
    std::vector<int> snp_num_in_gene(gene_num);
    std::map<std::string, int> snp_name_map;
    for (i = 0; i < snp_name.size(); i++) snp_name_map.insert(std::pair<std::string,int>(snp_name[i], i));
    for (i = 0; i < gene_num; i++) {
        auto iter1 = snp_name_map.find(gene2snp_1[i]);
        auto iter2 = snp_name_map.find(gene2snp_2[i]);
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
        for (j = iter1->second; j <= iter2->second; j++) chisq_o[i] += snp_chisq[j];
        Eigen::MatrixXf C;
        gbat_calcu_ld(X, sumsq_x, iter1->second, iter2->second, C); // iter2->second-1 because iter2 is one step further
        if(snp_num_in_gene[i] == 1) gene_pval[i] = StatFunc::pchisq(chisq_o[i], 1.0);
        else {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> saes(C);
            gene_pval[i] = StatFunc::pchisqsum(chisq_o[i], saes.eigenvalues().cast<double>());
        }

        if((i + 1) % 100 == 0 || (i + 1) == gene_num) LOGGER << i + 1 << " of " << gene_num << " genes.\r";
    }

    std::string filename = _out + ".gbat";
    LOGGER << "\nSaving the results of the gene-based association analysis to [" + filename + "] ..." << std::endl;
    std::ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Gene\tChr\tStart\tEnd\tNo.SNPs\tSNP_start\tSNP_end\tChisq(Obs)\tPvalue" << std::endl;
    for (i = 0; i < gene_num; i++) {
        if(gene_pval[i]>1.5) continue;
        ofile << gene_name[i] << "\t" << gene_chr[i] << "\t" << gene_bp1[i] << "\t" << gene_bp2[i] << "\t";
        ofile << snp_num_in_gene[i] << "\t" << gene2snp_1[i] << "\t" << gene2snp_2[i] << "\t" << chisq_o[i] << "\t" << gene_pval[i] << std::endl;
        //else ofile << "0\tNA\tNA\tNA\tNA" << std::endl;
    }
    ofile.close();
}

double gcta::gbat_simu_p(int &seed, int size, eigenMatrix &L, int simu_num, double chisq_o) {
    int i = 0, j = 0;
    std::vector<float> simu_chisq(simu_num);

    // debug
    LOGGER << "here simulation starts." << std::endl;

    //default_random_engine eng;
    // normal_distribution<float> rnorm(0.0, 1.0);


    #pragma omp parallel for private(j)
    for (i = 0; i < simu_num; i++) {
        eigenVector vec(size);
        //for(j=0; j<size; j++) vec[j]=rnorm(eng);

        // debug
        /*eigenVector tmp=L*vec;
        LOGGER<<"vec*L: "<<std::endl;
        LOGGER<<tmp<<std::endl;*/

        simu_chisq[i] = (L * vec).squaredNorm();
    }

    // debug
    /*LOGGER<<"chisq_o = "<<chisq_o<<std::endl;
    LOGGER<<"simu_chisq = ";
    for(i=0; i<simu_num; i++) LOGGER<<simu_chisq[i]<<" ";
    LOGGER<<std::endl;*/

    // debug
    LOGGER << "here find starts." << std::endl;


    int pos = (upper_bound(simu_chisq.begin(), simu_chisq.end(), chisq_o) - simu_chisq.begin());

    // debug
    LOGGER << "here find ends." << std::endl;


    // debug
    //LOGGER<<"pos = "<<pos<<std::endl;

    double pval = (double) (simu_num - pos) / (double) (simu_num);
    return (pval);
}
