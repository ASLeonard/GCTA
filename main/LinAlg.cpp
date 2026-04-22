/*
 *  GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Backend-agnostic linear algebra helpers (formerly MKL-specific).
 * All raw float/LAPACK code has been replaced with Eigen.
 *
 * 2014 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "gcta.h"

///////////
// reml

bool gcta::comput_inverse_logdet_LDLT(eigenMatrix &Vi, double &logdet) {
    Eigen::LLT<Eigen::MatrixXd> llt(Vi.template cast<double>());
    if (llt.info() != Eigen::Success) return false;
    logdet = 2.0 * llt.matrixL().toDenseMatrix().diagonal().array().log().sum();
    Vi = llt.solve(Eigen::MatrixXd::Identity(Vi.rows(), Vi.cols()))
             .template cast<typename eigenMatrix::Scalar>();
    return true;
}

bool gcta::comput_inverse_logdet_LU(eigenMatrix &Vi, double &logdet) {
    Eigen::FullPivLU<Eigen::MatrixXd> lu(Vi.template cast<double>());
    if (!lu.isInvertible()) return false;
    logdet = lu.matrixLU().diagonal().array().abs().log().sum();
    Vi = lu.inverse().template cast<typename eigenMatrix::Scalar>();
    return true;
}

bool gcta::comput_inverse_logdet_LU_array(int n, float *Vi, double &logdet) {
    Eigen::FullPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> lu(Eigen::Map<Eigen::Matrix<float,  Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(Vi, n, n).cast<double>());
    if (!lu.isInvertible()) return false;
    logdet = lu.matrixLU().diagonal().array().abs().log().sum();
    Eigen::Map<Eigen::Matrix<float,  Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(Vi, n, n) = lu.inverse().cast<float>();
    return true;
}

//////////////
// LD

// Compute r² between every pair of SNPs in each window block and flag for pruning.
// X must be already centred and divided by std (output of std_XMat with divid_by_std=true).
// r²(j,k) = (X[:,j]·X[:,k] / n)²  — symmetric matrix, so column-major Eigen storage
// is equivalent to row-major for rm_cor_snp which only reads the lower triangle.
void gcta::LD_pruning_blk(const Eigen::MatrixXf &X, std::vector<int> &brk_pnt, double rsq_cutoff, std::vector<int> &rm_snp_ID1)
{
    const long n = static_cast<long>(X.rows());

    for (size_t i = 0; i + 1 < brk_pnt.size(); i++) {
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        int size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        // n×size block; (X^T X)/n then squared → r² matrix (size×size)
        Eigen::MatrixXf X_sub = X.block(0, brk_pnt[i], n, size);
        Eigen::MatrixXf rsq = (X_sub.transpose() * X_sub) * (1.0f / n);
        rsq = rsq.array().square();

        std::vector<int> rm_snp_buf;
        rm_cor_snp(size, brk_pnt[i], rsq.data(), rsq_cutoff, rm_snp_buf);
        rm_snp_ID1.insert(rm_snp_ID1.end(), rm_snp_buf.begin(), rm_snp_buf.end());
    }
}

void gcta::LD_pruning(double rsq_cutoff, int wind_size) {
    check_autosome();

    make_XMat(_geno);
    eigenVector sd_SNP;
    std_XMat(_geno, sd_SNP, false, true, true);

    LOGGER << "\nPruning SNPs for LD ..." << std::endl;
    std::vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_size * 2);

    std::vector<int> rm_snp_indx;
    LD_pruning_blk(_geno, brk_pnt1, rsq_cutoff, rm_snp_indx);
    if (brk_pnt2.size() > 1) LD_pruning_blk(_geno, brk_pnt2, rsq_cutoff, rm_snp_indx);
    std::stable_sort(rm_snp_indx.begin(), rm_snp_indx.end());
    rm_snp_indx.erase(unique(rm_snp_indx.begin(), rm_snp_indx.end()), rm_snp_indx.end());

    int m_sub = static_cast<int>(rm_snp_indx.size());
    std::vector<std::string> rm_snp_name(m_sub);
    #pragma omp parallel for
    for (int idx = 0; idx < m_sub; idx++) rm_snp_name[idx] = _snp_name[_include[rm_snp_indx[idx]]];
    update_id_map_rm(rm_snp_name, _snp_name_map, _include);
    unsigned long m = _include.size();

    LOGGER << "After LD-pruning, " << m << " SNPs remain." << std::endl;
    std::string pruned_file = _out + ".prune.in";
    std::ofstream oprune(pruned_file.data());
    for (unsigned long i = 0; i < m; i++) oprune << _snp_name[_include[i]] << std::endl;
    oprune << std::endl;
    LOGGER << "The list of " << m << " LD-pruned SNPs (pruned in) has been saved in the file [" + pruned_file + "]." << std::endl;
}
