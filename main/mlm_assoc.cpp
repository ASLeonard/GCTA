/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for mixed linera model association analysis
 *
 * 2013 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "gcta.h"

#include <algorithm>
#include <fstream>
#include <string_view>

void gcta::mlma(std::string grm_file, bool m_grm_flag, std::string subtract_grm_file, std::string phen_file, std::string qcovar_file, std::string covar_file, int mphen, int MaxIter, std::vector<double> reml_priors, std::vector<double> reml_priors_var, bool no_constrain, bool within_family, bool inbred, bool no_adj_covar, std::string weight_file, std::string save_reml_file, std::string load_reml_file)
{
    _within_family=within_family;
    _reml_max_iter=MaxIter;
    unsigned long i = 0, j = 0, k = 0;
    bool grm_flag=(!grm_file.empty());
    bool qcovar_flag=(!qcovar_file.empty());
    bool covar_flag=(!covar_file.empty());
    if (!qcovar_flag && !covar_flag) no_adj_covar=false;
    if (m_grm_flag) grm_flag = false;
    bool subtract_grm_flag = (!subtract_grm_file.empty());
    const bool skip_grm_loading = !load_reml_file.empty();
    if (subtract_grm_flag && m_grm_flag) LOGGER.e(0, "the --mlma-subtract-grm option cannot be used in combination with the --mgrm option.");
    
    // Read data
    int qcovar_num=0, covar_num=0;
    std::vector<std::string> phen_ID, qcovar_ID, covar_ID, grm_id;
    std::vector< std::vector<std::string> > phen_buf, qcovar, covar; // save individuals by column
    std::vector<std::string> grm_files;
    
    if(phen_file.empty()){
        LOGGER.e(0, "no file name in --pheno.");
    }
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if(qcovar_flag){
        qcovar_num=read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if(covar_flag){
        covar_num=read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    // grm operations will overwrite the _keep
    if(_keep.size() < 1){
        LOGGER.e(0, "no individual is in common among the input files.");
    }

    if (skip_grm_loading) {
        if (!grm_file.empty() || m_grm_flag || subtract_grm_flag) {
            LOGGER << "Skipping GRM loading because --load-reml is set. Ensure the saved REML state matches the current sample set." << std::endl;
        }
    } else {
        if(subtract_grm_flag){
            grm_files.push_back(grm_file);
            grm_files.push_back(subtract_grm_file);
            for (i = 0; i < grm_files.size(); i++) {
                read_grm(grm_files[i], grm_id, false, true, true);
                update_id_map_kp(grm_id, _id_map, _keep);
            }
        }
        else{
            if(grm_flag){
                grm_files.push_back(grm_file);
                read_grm(grm_file, grm_id, true, false, true);
                update_id_map_kp(grm_id, _id_map, _keep);
            }
            else if (m_grm_flag) {
                read_grm_filenames(grm_file, grm_files, false);
                for (i = 0; i < grm_files.size(); i++) {
                    read_grm(grm_files[i], grm_id, false, true, true);
                    update_id_map_kp(grm_id, _id_map, _keep);
                }
            }
            else{
                grm_files.push_back("NA");
                make_grm(false, false, inbred, true, 0, true);
                for(i=0; i<_keep.size(); i++) grm_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
            }
        }
    }
    
    std::vector<std::string> uni_id;
    std::map<std::string, int> uni_id_map;
    make_uni_id(uni_id, uni_id_map);
    _n=_keep.size();
    if(_n<1) LOGGER.e(0, "no individual is in common in the input files.");
    LOGGER<<_n<<" individuals are in common in these files."<<std::endl;
    
    // construct model terms
    _y.setZero(_n);
    for(size_t i = 0; i < phen_ID.size(); ++i){
        if(auto iter = uni_id_map.find(phen_ID[i]); iter != uni_id_map.end()) {
            _y[iter->second] = std::stod(phen_buf[i][mphen-1]);
        }
    }

    _r_indx.clear();
    std::vector<int> kp;
    if (!skip_grm_loading) {
        if (subtract_grm_flag) {
            _r_indx = {0, 1};
            _A.resize(_r_indx.size());

            LOGGER << "\nReading the primary GRM from [" << grm_files[1] << "] ..." << std::endl;
            read_grm(grm_files[1], grm_id, true, false, false);

            StrFunc::match(uni_id, grm_id, kp);
            (_A[0]).resize(_n, _n);
            eigenMatrix A_N_buf(_n, _n);
            {
                eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
                Eigen::MatrixXf grm_N_sym_f(_grm_N.selfadjointView<Eigen::Lower>());
                eigenMatrix grm_N_sym = grm_N_sym_f.cast<double>();
                Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
                (_A[0]) = grm_sym(kp_idx, kp_idx);
                A_N_buf = grm_N_sym(kp_idx, kp_idx);
            }

            LOGGER << "\nReading the secondary GRM from [" << grm_files[0] << "] ..." << std::endl;
            read_grm(grm_files[0], grm_id, true, false, false);
            LOGGER<<"\nSubtracting [" << grm_files[1] << "] from [" << grm_files[0] << "] ..." << std::endl;
            StrFunc::match(uni_id, grm_id, kp);
            {
                eigenMatrix grm2_sym(_grm.selfadjointView<Eigen::Lower>());
                Eigen::MatrixXf grm2_N_sym_f(_grm_N.selfadjointView<Eigen::Lower>());
                eigenMatrix grm2_N_sym = grm2_N_sym_f.cast<double>();
                Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
                eigenMatrix grm2_slice = grm2_sym(kp_idx, kp_idx);
                eigenMatrix grm2_N_slice = grm2_N_sym(kp_idx, kp_idx);
                _A[0] = ((_A[0].array() * A_N_buf.array()) - (grm2_slice.array() * grm2_N_slice.array()))
                         / (A_N_buf.array() - grm2_N_slice.array());
            }
            _grm.resize(0,0);
            _grm_N.resize(0,0);
        }
        else {
            _r_indx.resize(grm_files.size() + 1);
            std::iota(_r_indx.begin(), _r_indx.end(), 0);
            _A.resize(_r_indx.size());
            if(grm_flag){
                StrFunc::match(uni_id, grm_id, kp);
                {
                    eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
                    Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
                    (_A[0]) = grm_sym(kp_idx, kp_idx);
                }
                _grm.resize(0,0);
            }
            else if(m_grm_flag){
                LOGGER << "There are " << grm_files.size() << " GRM file names specified in the file [" + grm_file + "]." << std::endl;
                for (i = 0; i < grm_files.size(); i++) {
                    LOGGER << "Reading the GRM from the " << i + 1 << "th file ..." << std::endl;
                    read_grm(grm_files[i], grm_id, true, false, true);
                    StrFunc::match(uni_id, grm_id, kp);
                    {
                        eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
                        Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
                        (_A[i]) = grm_sym(kp_idx, kp_idx);
                    }
                }
            }
            else{
                StrFunc::match(uni_id, grm_id, kp);
                {
                    eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
                    Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
                    (_A[0]) = grm_sym(kp_idx, kp_idx);
                }
                _grm.resize(0,0);
            }
        }
        _A[_r_indx.size()-1]=eigenMatrix::Identity(_n, _n);
        
        if(!weight_file.empty()){
            std::vector<std::string> weight_ID;
            std::vector<double> weights;

            read_weight(weight_file, weight_ID, weights);
            update_id_map_kp(weight_ID, _id_map, _keep);

        //if(!weight_file.empty()){
            // contruct weight
            Eigen::VectorXd v_weight(_n);
            for (size_t i = 0; i < weight_ID.size(); ++i) {
                if (auto it = uni_id_map.find(weight_ID[i]); it != uni_id_map.end()) {
                    v_weight(it->second) = weights[i];
                }
            }
            //v_weight = 1.0 / v_weight.array();
            //v_weight = 1.0 / v_weight.array() - (1.0 / v_weight.array()).mean() + 1;
            /*
            std::ofstream o_test("weight_out.txt");
            for(int i = 0; i < v_weight.size(); i++){
                o_test << v_weight[i] << "\t" << _y[i] << std::endl;
            }
            o_test.close();
            */

            _A[_r_indx.size() - 1].diagonal() = v_weight;
        }
    }

    // construct X matrix
    std::vector<eigenMatrix> E_float;
    eigenMatrix qE_float;
    construct_X(_n, uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);

    
    // names of variance component
    if (!skip_grm_loading) {
        for (size_t i = 0; i < grm_files.size(); ++i) {
            const std::string suffix = (grm_files.size() == 1) ? "" : std::to_string(i + 1);
            _var_name.emplace_back("V(G" + suffix + ")");
            _hsq_name.emplace_back("V(G" + suffix + ")/Vp");
        }
        _var_name.push_back("V(e)");
        
        // within family
        if(_within_family) detect_family();
    }
    
    // run REML algorithm
    LOGGER << "\nPerforming MLM association analyses" << (subtract_grm_flag?"":" (including the candidate SNP)") << " ..."<<std::endl;
    unsigned long n=_keep.size(), m=_include.size();
	
    if(!load_reml_file.empty()) {
            // Load REML state from file
            LOGGER << "Loading REML state from [" << load_reml_file << "] ..." << std::endl;
            load_reml_state(load_reml_file, no_adj_covar);
    } else {
        // Run REML estimation
        reml(false, true, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true);
        
        if(!save_reml_file.empty()) {
            // Save REML state and exit - no need for genotype data beyond this point
            LOGGER << "Saving REML state to [" << save_reml_file << "] ..." << std::endl;
            save_reml_state(save_reml_file, no_adj_covar);
            LOGGER << "REML estimation completed. Use --load-reml to perform association tests." << std::endl;
            
            // Clean up matrices that won't be needed
            _P.resize(0,0);
            _A.clear();
            return;
        }
    }
    
    _P.resize(0,0);
    _A.clear();
    LOGGER << "\nRegression coefficient(s) (e.g., general mean): " <<_b <<std::endl;

    std::vector<float> y(n);
    if(!no_adj_covar)
        Eigen::Map<Eigen::VectorXf>(y.data(), n) = (_y - _X*_b).cast<float>(); // adjust phenotype for covariates
    else
        Eigen::Map<Eigen::VectorXf>(y.data(), n) = _y.cast<float>();
    
/*    if(grm_flag || m_grm_flag){
        LOGGER<<std::endl;
        _geno_mkl=new float[n*m];
        make_XMat_mkl(_geno_mkl, false);
        #pragma omp parallel for private(j, k)
        for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                k=i*m+j;
                if(_geno_mkl[k]<1e5) _geno_mkl[k]-=_mu[_include[j]];
                else _geno_mkl[k]=0.0;
            }
        }
    }*/
    
    if (_mu.empty()) calcu_mu();

    //should we just init these here
    auto [beta, se, pval] = no_adj_covar
        ? mlma_calcu_stat_covar(std::span<const float>(y), m)
        : mlma_calcu_stat(std::span<const float>(y), m);
    _geno.resize(0, 0);
    
    const std::string filename = _out + ".mlma";
    LOGGER<<"\nSaving the results of the mixed linear model association analyses of "<<m<<" SNPs to ["+filename+"] ..."<<std::endl;
    std::ofstream ofile(filename);
    if(!ofile) LOGGER.e(0, "cannot open the file ["+filename+"] to write.");
    ofile<<"Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\t"<<(_log_pval ? "log_p" : "p")<<std::endl;
	for(size_t i = 0; i < m; ++i){
        const auto j = _include[i];
        ofile<<_chr[j]<<"\t"<<_snp_name[j]<<"\t"<<_bp[j]<<"\t"<<_ref_A[j]<<"\t"<<_other_A[j]<<"\t";
        if(pval[i]>1.5) ofile<<"NA\tNA\tNA\tNA"<<std::endl;
        else ofile<<0.5*_mu[j]<<"\t"<<beta[i]<<"\t"<<se[i]<<"\t"<<pval[i]<<std::endl;
    }
    ofile.close();
}

gcta::MlmaResult gcta::mlma_calcu_stat(std::span<const float> y, unsigned long m)
{
    const auto n = static_cast<Eigen::Index>(y.size());
    constexpr Eigen::Index max_block_size = 10000;
    unsigned long i = 0;
    std::vector<float> beta(m, 0.0f), se(m, 0.0f), pval(m, 2.0f);

    // Symmetrise Vi once into a plain matrix so every downstream GEMM/GEMV
    // dispatches through Eigen's full BLAS path rather than the slower
    // selfadjointView fallback that activates for block-expression targets.
    Eigen::MatrixXf Vi = _Vi.cast<float>().selfadjointView<Eigen::Upper>();
    _Vi.resize(0,0);

    // Precompute Vi * y once (plain GEMV, Vi is now fully dense symmetric)
    Eigen::Map<const Eigen::VectorXf> y_vec(y.data(), n);
    Eigen::VectorXf Vi_y(n);
    Vi_y.noalias() = Vi * y_vec;

    LOGGER<<"\nRunning association tests for "<<m<<" SNPs ..."<<std::endl;

    int k = 0, l = 0;
    Eigen::MatrixXf X_block, Vi_X_block(n, max_block_size);
    // Pre-allocated temporaries — no heap activity inside the hot loop
    Eigen::VectorXf Xt_Vi_y_block(max_block_size);
    Eigen::VectorXf xvx_diag(max_block_size);   // diag(X^T Vi X)
    std::vector<int> indx;

    int last_pct = -1;
    for(i = 0; i < m; ){
        int cur_pct = static_cast<int>(i * 100 / m);
        if(cur_pct != last_pct){
            LOGGER.p(0, std::to_string(i) + " / " + std::to_string(m) + " SNPs (" + std::to_string(cur_pct) + "%)");
            last_pct = cur_pct;
        }
        const Eigen::Index bs = static_cast<Eigen::Index>(
            std::min(static_cast<unsigned long>(max_block_size), m - i));
        indx.resize(bs);
        for(k = 0; k < bs; k++) indx[k] = i + k;
        make_XMat_subset(X_block, indx, false);  // X_block is n×bs

        // Vi * X_block → Vi_X_block (n×bs); Vi is plain dense → full GEMM
        Vi_X_block.leftCols(bs).noalias() = Vi * X_block;

        // diag(X^T Vi X): column-wise dot products via elementwise product + colwise sum
        // Avoids O(bs²) full matrix multiply; reads each column once.
        xvx_diag.head(bs) = (X_block.cwiseProduct(Vi_X_block.leftCols(bs))).colwise().sum();

        // X^T Vi_y (bs×1) — single GEMV, pre-allocated target
        Xt_Vi_y_block.head(bs).noalias() = X_block.transpose() * Vi_y;

        for(l = 0; l < bs; l++){
            const float inv_xvx = 1.0f / xvx_diag[l];
            if(std::isfinite(inv_xvx) && inv_xvx > 1.0e-30f){
                beta[i + l] = inv_xvx * Xt_Vi_y_block[l];
                se[i + l] = std::sqrt(inv_xvx);
                const float chisq = beta[i + l] / se[i + l];
                pval[i + l] = StatFunc::pchisq(chisq * chisq, 1, _log_pval);
            }
        }

        i += bs;
    }
    LOGGER.p(0, std::to_string(m) + " / " + std::to_string(m) + " SNPs (100%)");
    return {beta, se, pval};
}

gcta::MlmaResult gcta::mlma_calcu_stat_covar(std::span<const float> y, unsigned long m)
{
    const auto n = static_cast<Eigen::Index>(y.size());
    const Eigen::Index p = static_cast<Eigen::Index>(_X_c);  // number of fixed covariates
    constexpr Eigen::Index max_block_size = 10000;
    unsigned long i = 0;
    std::vector<float> beta(m, 0.0f), se(m, 0.0f), pval(m, 2.0f);

    // Symmetrise Vi once into a plain dense matrix so all downstream GEMMs
    // go through the full BLAS path (selfadjointView products on block-expression
    // targets may fall back to a scalar path in some Eigen/BLAS configurations).
    Eigen::MatrixXf Vi = _Vi.cast<float>().selfadjointView<Eigen::Upper>();
    _Vi.resize(0,0);

    // Precompute Vi * y once (plain GEMV on the now-dense Vi)
    Eigen::Map<const Eigen::VectorXf> y_vec(y.data(), n);
    Eigen::VectorXf Vi_y(n);
    Vi_y.noalias() = Vi * y_vec;

    // --- Schur complement precomputation (done once, not per SNP) ---
    // C (n×p): fixed covariate matrix
    Eigen::MatrixXf C = _X.cast<float>();

    // Vi_C = Vi * C  (n×p); plain GEMM, no block expressions
    Eigen::MatrixXf Vi_C(n, p);
    Vi_C.noalias() = Vi * C;

    // A = C^T Vi_C  (p×p).  Use Vi_C^T * C = C^T * Vi_C (same result, Vi symmetric).
    // Factor immediately; we never need the explicit inverse.
    Eigen::MatrixXf A_mat(p, p);
    A_mat.noalias() = Vi_C.transpose() * C;
    Vi_C.resize(0, 0);  // no longer needed — free n×p memory before the SNP loop
    Eigen::LLT<Eigen::MatrixXf> A_llt(A_mat);
    if(A_llt.info() != Eigen::Success)
        LOGGER.e(0, "covariate matrix C^T Vi C is not positive-definite/invertible.");

    // t = A^{-1} (C^T Vi_y)  (p×1)
    Eigen::VectorXf t_vec(p);
    t_vec.noalias() = C.transpose() * Vi_y;
    t_vec = A_llt.solve(t_vec);

    LOGGER << "\nRunning association tests for " << m << " SNPs ..." << std::endl;

    int k = 0, l = 0;
    Eigen::MatrixXf X_block;
    std::vector<int> indx;

    // All block temporaries pre-allocated outside the loop — zero heap activity per SNP
    Eigen::MatrixXf Vi_X_block(n, max_block_size);   // Vi * X_block
    Eigen::MatrixXf D_block(p, max_block_size);      // C^T * Vi_X_block  (p×bs)
    Eigen::MatrixXf E_block(p, max_block_size);      // A^{-1} * D_block  (p×bs)
    Eigen::VectorXf f_vec(max_block_size);           // X^T Vi_y          (bs×1)
    Eigen::VectorXf Dt_t_vec(max_block_size);        // D^T t             (bs×1)
    Eigen::VectorXf xvx_diag(max_block_size);        // diag(X^T Vi X)    (bs×1)
    Eigen::VectorXf d_dot_e_diag(max_block_size);    // diag(D^T E)       (bs×1)

    int last_pct_c = -1;
    for(i = 0; i < m; ){
        int cur_pct_c = static_cast<int>(i * 100 / m);
        if(cur_pct_c != last_pct_c){
            LOGGER.p(0, std::to_string(i) + " / " + std::to_string(m) + " SNPs (" + std::to_string(cur_pct_c) + "%)");
            last_pct_c = cur_pct_c;
        }
        const Eigen::Index bs = static_cast<Eigen::Index>(
            std::min(static_cast<unsigned long>(max_block_size), m - i));
        indx.resize(bs);
        for(k = 0; k < bs; k++) indx[k] = i + k;
        make_XMat_subset(X_block, indx, false);  // X_block is n×bs

        // Vi * X_block → Vi_X_block (n×bs); plain GEMM
        Vi_X_block.leftCols(bs).noalias() = Vi * X_block;

        // D = C^T * Vi_X_block  (p×bs); plain GEMM on evaluated sub-matrices
        D_block.leftCols(bs).noalias() = C.transpose() * Vi_X_block.leftCols(bs);

        // E = A^{-1} * D  (p×bs) via LLT; noalias safe since D and E are distinct
        E_block.leftCols(bs).noalias() = A_llt.solve(D_block.leftCols(bs));

        // f = X^T Vi_y  (bs×1)
        f_vec.head(bs).noalias() = X_block.transpose() * Vi_y;

        // Dt_t = D^T t  (bs×1)
        Dt_t_vec.head(bs).noalias() = D_block.leftCols(bs).transpose() * t_vec;

        // diag(X^T Vi X) and diag(D^T E) via elementwise product + colwise sum;
        // avoids forming the full bs×bs products.
        xvx_diag.head(bs)    = (X_block.cwiseProduct(Vi_X_block.leftCols(bs))).colwise().sum();
        d_dot_e_diag.head(bs) = (D_block.leftCols(bs).cwiseProduct(E_block.leftCols(bs))).colwise().sum();

        for(l = 0; l < bs; l++){
            const float S = xvx_diag[l] - d_dot_e_diag[l];  // Schur complement = 1/Var(beta_snp)
            beta[i+l] = (f_vec[l] - Dt_t_vec[l]) / S;
            if(S > 1.0e-30f){
                se[i+l] = std::sqrt(1.0f / S);
                const float chisq = beta[i+l] / se[i+l];
                pval[i+l] = StatFunc::pchisq(chisq * chisq, 1, _log_pval);
            }
        }

        i += bs;
    }
    LOGGER.p(0, std::to_string(m) + " / " + std::to_string(m) + " SNPs (100%)");
    return {beta, se, pval};
}

void gcta::mlma_loco(std::string phen_file, std::string qcovar_file, std::string covar_file, int mphen, int MaxIter, std::vector<double> reml_priors, std::vector<double> reml_priors_var, bool no_constrain, bool inbred, bool no_adj_covar)
{
    unsigned long i=0, j=0, k=0, c1=0, c2=0, n=0;
    _reml_max_iter=MaxIter;
    bool qcovar_flag=(!qcovar_file.empty());
    bool covar_flag=(!covar_file.empty());
    if(!qcovar_flag && !covar_flag) no_adj_covar=false;
    
    // Read data
    int qcovar_num=0, covar_num=0;
    std::vector<std::string> phen_ID, qcovar_ID, covar_ID, grm_id;
    std::vector< std::vector<std::string> > phen_buf, qcovar, covar; // save individuals by column

    if(phen_file.empty()){
        LOGGER.e(0, "no file name in --pheno.");
    }
    
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if(qcovar_flag){
        qcovar_num=read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if(covar_flag){
        covar_num=read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    n=_keep.size();
    _n=_keep.size();
    if(_n<1) LOGGER.e(0, "no individual is in common among the input files.");
    LOGGER<<_n<<" individuals are in common in these files."<<std::endl;
    
    std::vector<int> chrs, vi_buf(_chr);
    std::ranges::sort(vi_buf);
    auto unique_range = std::ranges::unique(vi_buf);
    vi_buf.erase(unique_range.begin(), vi_buf.end());
    if(vi_buf.size()<2) LOGGER.e(0, "There is only one chromosome. The MLM leave-on-chromosome-out (LOCO) analysis requires at least two chromosomes.");
    for(i=0; i<vi_buf.size(); i++){
        if(vi_buf[i]<=_autosome_num) chrs.push_back(vi_buf[i]);
    }
    std::vector<int> include_o(_include);
    std::map<std::string, int> snp_name_map_o(_snp_name_map);
    std::vector<float> m_chrs_f(chrs.size());
    std::vector<eigenMatrix> grm_chrs(chrs.size());
    std::vector< std::vector<int> > icld_chrs(chrs.size());
    LOGGER<<std::endl;
    if(_mu.empty()) calcu_mu();
    LOGGER<<"\nCalculating the genetic relationship matrix for each of the "<<chrs.size()<<" chromosomes ... "<<std::endl;
    for(c1=0; c1<chrs.size(); c1++){
        LOGGER<<"Chr "<<chrs[c1]<<":"<<std::endl;
        extract_chr(chrs[c1], chrs[c1]);
        make_grm(false, false, inbred, true, 0, true);
        
        m_chrs_f[c1]=(float)_include.size();
        icld_chrs[c1]=_include;
        _include=include_o;
        _snp_name_map=snp_name_map_o;
        
        grm_chrs[c1]=std::move(_grm);
        _geno.resize(0, 0);
    }
    for(i=0; i<_keep.size(); i++) grm_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
    
    std::vector<std::string> uni_id;
    std::map<std::string, int> uni_id_map;
    std::map<std::string, int>::iterator iter;
    make_uni_id(uni_id, uni_id_map);
    
    // construct model terms
    _y.setZero(_n);
    for(i=0; i<phen_ID.size(); i++){
        iter=uni_id_map.find(phen_ID[i]);
        if(iter==uni_id_map.end()) continue;
        _y[iter->second]=atof(phen_buf[i][mphen-1].c_str());
    }
    
    // construct X matrix
    std::vector<eigenMatrix> E_float;
    eigenMatrix qE_float;
    construct_X(_n, uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);
    
    // names of variance component
    _var_name.push_back("V(G)");
    _hsq_name.push_back("V(G)/Vp");
    _var_name.push_back("V(e)");
    
    // MLM association
    LOGGER<<"\nPerforming MLM association analyses (leave-one-chromosome-out) ..."<<std::endl;
    
    std::vector<int> kp;
    StrFunc::match(uni_id, grm_id, kp);
    _r_indx.resize(2);
    for(i=0; i<2; i++) _r_indx[i]=i;
    _A.resize(_r_indx.size());
    _A[1]=eigenMatrix::Identity(_n, _n);
    
    eigenVector y_buf=_y;
    std::vector<float> y(_n);
    std::vector<std::vector<float>> beta(chrs.size()), se(chrs.size()), pval(chrs.size());
    for(c1=0; c1<chrs.size(); c1++){
        LOGGER<<"\n-----------------------------------\n#Chr "<<chrs[c1]<<":"<<std::endl;
        extract_chr(chrs[c1], chrs[c1]);
        
        _A[0] = eigenMatrix::Zero(_n, _n);
        double d_buf = 0;
        {
            Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
            for(c2=0; c2<chrs.size(); c2++){
                if(chrs[c1]==chrs[c2]) continue;
                eigenMatrix grm_sym(grm_chrs[c2].selfadjointView<Eigen::Lower>());
                _A[0] += grm_sym(kp_idx, kp_idx) * static_cast<double>(m_chrs_f[c2]);
                d_buf += m_chrs_f[c2];
            }
        }
        _A[0] /= d_buf;
        
        // run REML algorithm
        reml(false, true, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true);
        if(!no_adj_covar) y_buf=_y.array()-(_X*_b).array(); // adjust phenotype for covariates
        for(i=0; i<_n; i++) y[i]=y_buf[i];
        reml_priors.clear();
        reml_priors_var=_varcmp;
        _P.resize(0,0);
        _A[0].resize(0,0);

        auto [b, s, p] = no_adj_covar
            ? mlma_calcu_stat_covar(std::span<const float>(y), _include.size())
            : mlma_calcu_stat(std::span<const float>(y), _include.size());
        beta[c1] = std::move(b); se[c1] = std::move(s); pval[c1] = std::move(p);
        
        _include=include_o;
        _snp_name_map=snp_name_map_o;
        LOGGER<<"-----------------------------------"<<std::endl;
    }
    
    std::string filename=_out+".loco.mlma";
    LOGGER<<"\nSaving the results of the mixed linear model association analyses of "<<_include.size()<<" SNPs to ["+filename+"] ..."<<std::endl;
    std::ofstream ofile(filename.c_str());
    if(!ofile) LOGGER.e(0, "cannot open the file ["+filename+"] to write.");
    ofile<<"Chr\tSNP\tbp\tA1\tA2\tFreq\tb\tse\t"<<(_log_pval ? "log_p" : "p")<<std::endl;
    for(c1=0; c1<chrs.size(); c1++){
        for(i=0; i<icld_chrs[c1].size(); i++){
            j=icld_chrs[c1][i];
            ofile<<_chr[j]<<"\t"<<_snp_name[j]<<"\t"<<_bp[j]<<"\t"<<_ref_A[j]<<"\t"<<_other_A[j]<<"\t";
            if(pval[c1][i]>1.5) ofile<<"NA\tNA\tNA\tNA"<<std::endl;
            else ofile<<0.5*_mu[j]<<"\t"<<beta[c1][i]<<"\t"<<se[c1][i]<<"\t"<<pval[c1][i]<<std::endl;
        }
    }
    ofile.close();
}

void gcta::grm_minus_grm(float *grm, float *sub_grm)
{
    int i=0, j=0, k=0, n=_n;
    
    #pragma omp parallel for private(j,k)
    for(i=0; i<n; i++){
		for(j=0; j<=i; j++){
            k=i*n+j;
            sub_grm[k]=grm[k]-sub_grm[k];
		}
	}

}

void gcta::save_reml_state(const std::string& filename, bool no_adj_covar)
{
    // Ensure _Vi is materialised regardless of which REML path was taken:
    //   - Default path:  _Vi_use_llt == false, _Vi is already an n×n matrix.
    //   - Hutch++ path:  _Vi_use_llt == true, _Vi was released to save RAM.
    //     _Vi_llt holds L from V = L Lᵀ so _Vi_llt.solve(I) gives V^{-1} in O(n²).
    if (_Vi_use_llt) {
        _Vi = _Vi_llt.solve(eigenMatrix::Identity(_n, _n));
        _Vi_use_llt = false;
    }

    std::ofstream outfile(filename, std::ios::binary);
    if(!outfile.is_open()) LOGGER.e(0, "cannot open the file ["+filename+"] to write.");

    // ---- Header ----
    // Magic "GOBY": packed-triangle symmetric format (replaces old "REML" dense format).
    struct Header {
        char    magic[4]    = {'G','O','B','Y'};
        int32_t n           = 0;
        int32_t x_c         = 0;
        int32_t num_varcmp  = 0;
        int32_t num_r_indx  = 0;
    } hdr;
    hdr.n          = static_cast<int32_t>(_n);
    hdr.x_c        = static_cast<int32_t>(_X_c);
    hdr.num_varcmp = static_cast<int32_t>(_varcmp.size());
    hdr.num_r_indx = static_cast<int32_t>(_r_indx.size());
    outfile.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    // ---- _Vi: packed lower triangle (float) ----
    // _Vi is symmetric positive-definite, so only the n(n+1)/2 lower-triangle
    // elements are unique.  We save those directly rather than the full n×n matrix.
    // This halves the file size vs the old dense layout.
    //
    // We save _Vi (the inverse) rather than the LLT factor L for two reasons:
    //   1. On the default (non-Hutch++) path _Vi_llt is not reliably populated
    //      after convergence — only _Vi is guaranteed to exist.
    //   2. Serialising L as float and reconstructing L Lᵀ introduces squared
    //      rounding errors before the re-factorisation, degrading numerical quality.
    {
        const int32_t n = hdr.n;
        const size_t tri_size = static_cast<size_t>(n) * (n + 1) / 2;
        std::vector<float> tri(tri_size);
        size_t idx = 0;
        for(int32_t j = 0; j < n; ++j)
            for(int32_t i = j; i < n; ++i)
                tri[idx++] = static_cast<float>(_Vi(i, j));
        outfile.write(reinterpret_cast<const char*>(tri.data()),
                      static_cast<std::streamsize>(tri_size * sizeof(float)));
    }

    // ---- _b (float, only when covariates are present) ----
    if(!no_adj_covar) {
        Eigen::VectorXf b_f = _b.cast<float>();
        outfile.write(reinterpret_cast<const char*>(b_f.data()),
                      static_cast<std::streamsize>(hdr.x_c * sizeof(float)));
    }

    // ---- _varcmp (float) ----
    {
        Eigen::VectorXf vc_f =
            Eigen::Map<const Eigen::VectorXd>(_varcmp.data(), hdr.num_varcmp).cast<float>();
        outfile.write(reinterpret_cast<const char*>(vc_f.data()),
                      static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));
    }

    // ---- _r_indx (int32) ----
    outfile.write(reinterpret_cast<const char*>(_r_indx.data()),
                  static_cast<std::streamsize>(hdr.num_r_indx * sizeof(int)));

    outfile.flush();
    if(!outfile)
        LOGGER.e(0, "write error on ["+filename+"]: disk full or I/O failure.");

    const size_t tri_bytes  = static_cast<size_t>(hdr.n) * (hdr.n + 1) / 2 * sizeof(float);
    const size_t full_bytes = static_cast<size_t>(hdr.n) * hdr.n * sizeof(float);
    LOGGER << "Saved REML state (n=" << _n << ", covariates=" << _X_c
           << ", variance components=" << hdr.num_varcmp << ") to [" << filename << "] ("
           << tri_bytes / (1<<20) << " MiB vs " << full_bytes / (1<<20) << " MiB dense)." << std::endl;
}

void gcta::load_reml_state(const std::string& filename, bool no_adj_covar)
{
    std::ifstream infile(filename, std::ios::binary);
    if(!infile.is_open()) LOGGER.e(0, "cannot open the file ["+filename+"] to read. "
        "Make sure you have run --mlma --save-reml first with matching --out prefix.");

    auto must_read = [&](void* dst, std::streamsize n){
        infile.read(reinterpret_cast<char*>(dst), n);
        if(infile.gcount() != n)
            LOGGER.e(0, "unexpected end of file reading ["+filename+"] — file may be truncated.");
    };

    // ---- Header ----
    struct Header {
        char    magic[4];
        int32_t n, x_c, num_varcmp, num_r_indx;
    } hdr;
    must_read(&hdr, sizeof(hdr));

    if(std::string_view(hdr.magic, 4) != "GOBY")
        LOGGER.e(0, "file ["+filename+"] is not a valid REML state file.");
    if(hdr.n <= 0 || hdr.x_c < 0 || hdr.num_varcmp <= 0 || hdr.num_r_indx <= 0)
        LOGGER.e(0, "file ["+filename+"] contains invalid dimensions — file may be corrupt.");
    if(hdr.n != static_cast<int32_t>(_n))
        LOGGER.e(0, "sample size mismatch: REML state has n=" + std::to_string(hdr.n)
                  + " but current dataset has n=" + std::to_string(_n));
    if(hdr.x_c != static_cast<int32_t>(_X_c))
        LOGGER.e(0, "number of covariates mismatch: REML state has "
                  + std::to_string(hdr.x_c) + " but current dataset has " + std::to_string(_X_c));

    // ---- _Vi: packed lower triangle → full symmetric matrix ----
    // Read float triangle, upcast to double, and reconstruct the full n×n _Vi
    // by reflecting the lower triangle into the upper half.
    {
        const int32_t n = hdr.n;
        const size_t tri_size = static_cast<size_t>(n) * (n + 1) / 2;
        std::vector<float> tri(tri_size);
        must_read(tri.data(), static_cast<std::streamsize>(tri_size * sizeof(float)));

        _Vi.resize(n, n);
        size_t idx = 0;
        for(int32_t j = 0; j < n; ++j) {
            for(int32_t i = j; i < n; ++i) {
                const double v = static_cast<double>(tri[idx++]);
                _Vi(i, j) = v;
                _Vi(j, i) = v;   // reflect to upper triangle
            }
        }
        _Vi_use_llt = false;
    }

    // ---- _b ----
    if(!no_adj_covar) {
        Eigen::VectorXf b_f(hdr.x_c);
        must_read(b_f.data(), static_cast<std::streamsize>(hdr.x_c * sizeof(float)));
        _b = b_f.cast<double>();
    }

    // ---- _varcmp ----
    {
        Eigen::VectorXf vc_f(hdr.num_varcmp);
        must_read(vc_f.data(), static_cast<std::streamsize>(hdr.num_varcmp * sizeof(float)));
        _varcmp.resize(hdr.num_varcmp);
        Eigen::Map<Eigen::VectorXd>(_varcmp.data(), hdr.num_varcmp) = vc_f.cast<double>();
    }

    // ---- _r_indx ----
    _r_indx.resize(hdr.num_r_indx);
    must_read(_r_indx.data(),
              static_cast<std::streamsize>(hdr.num_r_indx * sizeof(int)));

    LOGGER << "Loaded REML state from [" << filename << "]: n=" << hdr.n
           << ", covariates=" << hdr.x_c
           << ", variance components=" << hdr.num_varcmp << std::endl;
}