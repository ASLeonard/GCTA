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
            LOGGER << "Skipping GRM loading because --load-reml is std::set. Ensure the saved REML state matches the current sample std::set." << std::endl;
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
                make_grm_mkl(false, false, inbred, true, 0, true);
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
            Eigen::MatrixXf A_N_buf(_n, _n);
            #pragma omp parallel for private(k)
            for (j = 0; j < _n; j++) {
                for (k = 0; k <= j; k++) {
                    if (kp[j] >= kp[k]){
                        (_A[0])(k, j) = (_A[0])(j, k) = _grm(kp[j], kp[k]);
                        A_N_buf(k, j) = A_N_buf(j, k) = _grm_N(kp[j], kp[k]);
                    }
                    else{
                        (_A[0])(k, j) = (_A[0])(j, k) = _grm(kp[k], kp[j]);
                        A_N_buf(k, j) = A_N_buf(j, k) = _grm_N(kp[k], kp[j]);
                    }
                }
            }

            LOGGER << "\nReading the secondary GRM from [" << grm_files[0] << "] ..." << std::endl;
            read_grm(grm_files[0], grm_id, true, false, false);
            LOGGER<<"\nSubtracting [" << grm_files[1] << "] from [" << grm_files[0] << "] ..." << std::endl;
            StrFunc::match(uni_id, grm_id, kp);
            #pragma omp parallel for private(k)
            for (j = 0; j < _n; j++) {
                for (k = 0; k <= j; k++) {
                    if (kp[j] >= kp[k]) (_A[0])(k, j) = (_A[0])(j, k) = ((_A[0])(j, k) * A_N_buf(j, k)  - _grm(kp[j], kp[k]) * _grm_N(kp[j], kp[k])) / (A_N_buf(j, k) - _grm_N(kp[j], kp[k]));
                    else (_A[0])(k, j) = (_A[0])(j, k) = ((_A[0])(j, k) * A_N_buf(j, k) - _grm(kp[k], kp[j]) * _grm_N(kp[k], kp[j])) / (A_N_buf(j, k) - _grm_N(kp[k], kp[j]));
                }
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
                (_A[0]).resize(_n, _n);
                #pragma omp parallel for private(j)
                for(i=0; i<_n; i++){
                    for(j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=_grm(kp[i],kp[j]);
                }
                _grm.resize(0,0);
            }
            else if(m_grm_flag){
                LOGGER << "There are " << grm_files.size() << " GRM file names specified in the file [" + grm_file + "]." << std::endl;
                for (i = 0; i < grm_files.size(); i++) {
                    LOGGER << "Reading the GRM from the " << i + 1 << "th file ..." << std::endl;
                    read_grm(grm_files[i], grm_id, true, false, true);
                    StrFunc::match(uni_id, grm_id, kp);
                    (_A[i]).resize(_n, _n);
                    #pragma omp parallel for private(k)
                    for (j = 0; j < _n; j++) {
                        for (k = 0; k <= j; k++) {
                            if (kp[j] >= kp[k]) (_A[i])(k, j) = (_A[i])(j, k) = _grm(kp[j], kp[k]);
                            else (_A[i])(k, j) = (_A[i])(j, k) = _grm(kp[k], kp[j]);
                        }
                    }
                }
            }
            else{
                StrFunc::match(uni_id, grm_id, kp);
                (_A[0]).resize(_n, _n);
                #pragma omp parallel for private(j)
                for(i=0; i<_n; i++){
                    for(j=0; j<=i; j++) (_A[0])(j,i)=(_A[0])(i,j)=_grm_mkl[kp[i]*_n+kp[j]];
                }
                delete[] _grm_mkl;
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
    std::cout << "\nRegression coefficient(s) (e.g., general mean): " <<_b <<std::endl;

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
    eigenVector beta, se, pval;
    if(no_adj_covar) mlma_calcu_stat_covar(std::span<const float>(y), std::span<const float>(_geno_mkl, static_cast<size_t>(n) * m), m, beta, se, pval);
    else mlma_calcu_stat(std::span<const float>(y), std::span<const float>(_geno_mkl, static_cast<size_t>(n) * m), m, beta, se, pval);
    delete[] _geno_mkl;
    
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

void gcta::mlma_calcu_stat(std::span<const float> y, [[maybe_unused]] std::span<const float> geno_mkl, unsigned long m, eigenVector &beta, eigenVector &se, eigenVector &pval)
{
    const auto n = static_cast<unsigned long>(y.size());
    constexpr int max_block_size = 10000;
    unsigned long i = 0, j = 0;
    std::vector<float> Vi(static_cast<size_t>(n) * n);
    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        for(j=0; j<n; j++) Vi[i*n+j]=_Vi(i,j);
    }
    _Vi.resize(0,0);

    // Precompute Vi * y once (Vi is symmetric → use ssymv)
    std::vector<float> Vi_y(n);
    cblas_ssymv(CblasColMajor, CblasUpper,
                n, 1.0f, Vi.data(), n,
                y.data(), 1,
                0.0f, Vi_y.data(), 1);

    beta.resize(m);
    se=eigenVector::Zero(m);
    pval=eigenVector::Constant(m,2);
    LOGGER<<"\nRunning association tests for "<<m<<" SNPs ..."<<std::endl;

    int k = 0, l = 0;
    Eigen::MatrixXf X_block;
    std::vector<int> indx;
    // Vi_X_block: col-major n × bs result of Vi * X_block
    std::vector<float> Vi_X_block(static_cast<size_t>(n) * max_block_size);
    std::vector<float> Xt_Vi_X_block(max_block_size);
    std::vector<float> Xt_Vi_y_block(max_block_size);

    for(i = 0; i < m; ){
        int bs = (int)std::min((unsigned long)max_block_size, m - i);
        indx.resize(bs);
        for(k = 0; k < bs; k++) indx[k] = i + k;
        make_XMat_subset(X_block, indx, false);

        // Vi (symmetric n×n) * X_block (col-major n×bs) → Vi_X_block (col-major n×bs)
        cblas_ssymm(CblasColMajor, CblasLeft, CblasUpper,
                    n, bs, 1.0f, Vi.data(), n,
                    X_block.data(), n,
                    0.0f, Vi_X_block.data(), n);

        // Diagonal of X_block^T * Vi_X_block: contiguous stride-1 dot products
        for(l = 0; l < bs; l++){
            Xt_Vi_X_block[l] = cblas_sdot(n, X_block.data() + (size_t)l * n, 1,
                                           Vi_X_block.data() + (size_t)l * n, 1);
        }

        // X_block^T * Vi_y (precomputed) — avoids re-reading n×bs Vi_X_block
        cblas_sgemv(CblasColMajor, CblasTrans,
                    n, bs,
                    1.0f, X_block.data(), n,
                    Vi_y.data(), 1,
                    0.0f, Xt_Vi_y_block.data(), 1);

        for(l = 0; l < bs; l++){
            float inv_xvx = 1.0f / Xt_Vi_X_block[l];
            se[i + l] = inv_xvx;
            beta[i + l] = inv_xvx * Xt_Vi_y_block[l];
            if(inv_xvx > 1.0e-30f){
                se[i + l] = std::sqrt(inv_xvx);
                float chisq = beta[i + l] / se[i + l];
                pval[i + l] = StatFunc::pchisq(chisq * chisq, 1, _log_pval);
            }
        }

        i += bs;
    }
}

void gcta::mlma_calcu_stat_covar(std::span<const float> y, [[maybe_unused]] std::span<const float> geno_mkl, unsigned long m, eigenVector &beta, eigenVector &se, eigenVector &pval)
{
    const auto n = static_cast<unsigned long>(y.size());
    const int p = static_cast<int>(_X_c);  // number of fixed covariates
    constexpr int max_block_size = 10000;
    unsigned long i = 0, j = 0;

    // Flatten Vi to a float array (Vi is symmetric so row-major == col-major)
    std::vector<float> Vi(static_cast<size_t>(n) * n);
    #pragma omp parallel for private(j)
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++) Vi[i*n+j] = static_cast<float>(_Vi(i,j));
    _Vi.resize(0,0);

    // Precompute Vi * y once (Vi is symmetric → use ssymv)
    std::vector<float> Vi_y(n);
    cblas_ssymv(CblasColMajor, CblasUpper,
                n, 1.0f, Vi.data(), n,
                y.data(), 1, 0.0f, Vi_y.data(), 1);

    // --- Schur complement precomputation (done once, not per SNP) ---
    // C (n×p, col-major float): fixed covariate matrix
    std::vector<float> C(static_cast<size_t>(n) * p);
    for(i = 0; i < n; i++)
        for(j = 0; j < static_cast<unsigned long>(p); j++)
            C[j*n + i] = static_cast<float>(_X(i,j));

    // Vi_C = Vi * C  (n×p, col-major; Vi is symmetric → use ssymm)
    std::vector<float> Vi_C(static_cast<size_t>(n) * p);
    cblas_ssymm(CblasColMajor, CblasLeft, CblasUpper,
                n, p, 1.0f, Vi.data(), n, C.data(), n, 0.0f, Vi_C.data(), n);

    // A = C^T Vi_C  (p×p), then invert in-place → A becomes A^{-1}
    std::vector<float> A(static_cast<size_t>(p) * p);
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                p, p, n, 1.0f, C.data(), n, Vi_C.data(), n, 0.0f, A.data(), p);
    double d_buf = 0.0;
    if(!comput_inverse_logdet_LU_mkl_array(p, A.data(), d_buf))
        LOGGER.e(0, "covariate matrix C^T Vi C is not invertible.");

    // c = C^T Vi_y  (p×1)
    std::vector<float> c_vec(p);
    cblas_sgemv(CblasColMajor, CblasTrans,
                n, p, 1.0f, C.data(), n, Vi_y.data(), 1, 0.0f, c_vec.data(), 1);

    // t = A^{-1} c  (p×1): used to adjust numerator of each SNP test
    std::vector<float> t_vec(p);
    cblas_sgemv(CblasColMajor, CblasNoTrans,
                p, p, 1.0f, A.data(), p, c_vec.data(), 1, 0.0f, t_vec.data(), 1);

    beta.resize(m);
    se = eigenVector::Zero(m);
    pval = eigenVector::Constant(m, 2);
    LOGGER << "\nRunning association tests for " << m << " SNPs ..." << std::endl;

    int k = 0, l = 0;
    Eigen::MatrixXf X_block;
    std::vector<int> indx;
    // Pre-allocate block buffers (sized for max block; reused each iteration)
    std::vector<float> Vi_X_block(static_cast<size_t>(n) * max_block_size);
    std::vector<float> D_block(static_cast<size_t>(p) * max_block_size);  // C^T * Vi_X_block (p×bs)
    std::vector<float> E_block(static_cast<size_t>(p) * max_block_size);  // A^{-1} * D_block
    std::vector<float> f_vec(max_block_size);     // X_block^T Vi_y
    std::vector<float> Dt_t_vec(max_block_size);  // D_block^T t

    for(i = 0; i < m; ){
        int bs = static_cast<int>(std::min(static_cast<unsigned long>(max_block_size), m - i));
        indx.resize(bs);
        for(k = 0; k < bs; k++) indx[k] = i + k;
        make_XMat_subset(X_block, indx, false);

        // Vi * X_block → Vi_X_block  (n×bs, col-major; Vi symmetric → ssymm)
        cblas_ssymm(CblasColMajor, CblasLeft, CblasUpper,
                    n, bs, 1.0f, Vi.data(), n, X_block.data(), n, 0.0f, Vi_X_block.data(), n);

        // D = C^T * Vi_X_block  (p×bs) — reuses Vi_X_block, avoids a second read of Vi_C
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    p, bs, n, 1.0f, C.data(), n, Vi_X_block.data(), n, 0.0f, D_block.data(), p);

        // E = A^{-1} * D  (p×bs)
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    p, bs, p, 1.0f, A.data(), p, D_block.data(), p, 0.0f, E_block.data(), p);

        // f = X_block^T Vi_y  (bs×1)
        cblas_sgemv(CblasColMajor, CblasTrans,
                    n, bs, 1.0f, X_block.data(), n, Vi_y.data(), 1, 0.0f, f_vec.data(), 1);

        // Dt_t = D_block^T t  (bs×1): covariate correction to numerator
        cblas_sgemv(CblasColMajor, CblasTrans,
                    p, bs, 1.0f, D_block.data(), p, t_vec.data(), 1, 0.0f, Dt_t_vec.data(), 1);

        for(l = 0; l < bs; l++){
            // a = x^T Vi x  (diagonal element of X_block^T Vi_X_block)
            float a = cblas_sdot(n, X_block.data()   + static_cast<size_t>(l)*n, 1,
                                    Vi_X_block.data() + static_cast<size_t>(l)*n, 1);
            // d^T e = D[:,l]^T E[:,l]  (covariate correction to denominator)
            float d_dot_e = cblas_sdot(p, D_block.data() + static_cast<size_t>(l)*p, 1,
                                          E_block.data() + static_cast<size_t>(l)*p, 1);
            float S = a - d_dot_e;  // Schur complement = 1/Var(beta_snp)
            beta[i+l] = (f_vec[l] - Dt_t_vec[l]) / S;
            if(S > 1.0e-30f){
                se[i+l] = std::sqrt(1.0f / S);
                float chisq = beta[i+l] / se[i+l];
                pval[i+l] = StatFunc::pchisq(chisq * chisq, 1, _log_pval);
            }
        }

        i += bs;
    }
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
    std::vector<float *> grm_chrs(chrs.size());
    std::vector<float *> geno_chrs(chrs.size());
    std::vector< std::vector<int> > icld_chrs(chrs.size());
    LOGGER<<std::endl;
    if(_mu.empty()) calcu_mu();
    LOGGER<<"\nCalculating the genetic relationship matrix for each of the "<<chrs.size()<<" chromosomes ... "<<std::endl;
    for(c1=0; c1<chrs.size(); c1++){
        LOGGER<<"Chr "<<chrs[c1]<<":"<<std::endl;
        extract_chr(chrs[c1], chrs[c1]);
        make_grm_mkl(false, false, inbred, true, 0, true);
        
        m_chrs_f[c1]=(float)_include.size();
        icld_chrs[c1]=_include;
        _include=include_o;
        _snp_name_map=snp_name_map_o;
        
        geno_chrs[c1]=_geno_mkl;
        _geno_mkl=NULL;
        grm_chrs[c1]=_grm_mkl;
        _grm_mkl=NULL;
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
    std::vector<eigenVector> beta(chrs.size()), se(chrs.size()), pval(chrs.size());
    for(c1=0; c1<chrs.size(); c1++){
        LOGGER<<"\n-----------------------------------\n#Chr "<<chrs[c1]<<":"<<std::endl;
        extract_chr(chrs[c1], chrs[c1]);
        
        _A[0]=eigenMatrix::Zero(_n, _n);
        double d_buf=0;
        for(c2=0; c2<chrs.size(); c2++){
            if(chrs[c1]==chrs[c2]) continue;
            #pragma omp parallel for private(j)
            for(i=0; i<_n; i++){
                for(j=0; j<=i; j++){
                    (_A[0])(i,j)+=(grm_chrs[c2])[kp[i]*_n+kp[j]]*m_chrs_f[c2];
                }
            }
            d_buf+=m_chrs_f[c2];
        }
        
        #pragma omp parallel for private(j)
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++){
                (_A[0])(i,j)/=d_buf;
                (_A[0])(j,i)=(_A[0])(i,j);
            }
        }
        
        // run REML algorithm
        reml(false, true, true, reml_priors, reml_priors_var, -2.0, -2.0, no_constrain, true, true);
        if(!no_adj_covar) y_buf=_y.array()-(_X*_b).array(); // adjust phenotype for covariates
        for(i=0; i<_n; i++) y[i]=y_buf[i];
        reml_priors.clear();
        reml_priors_var=_varcmp;
        _P.resize(0,0);
        _A[0].resize(0,0);

        if(no_adj_covar)  mlma_calcu_stat_covar(std::span<const float>(y), std::span<const float>(geno_chrs[c1], static_cast<size_t>(n) * _include.size()), _include.size(), beta[c1], se[c1], pval[c1]);
        else mlma_calcu_stat(std::span<const float>(y), std::span<const float>(geno_chrs[c1], static_cast<size_t>(n) * _include.size()), _include.size(), beta[c1], se[c1], pval[c1]);
        
        _include=include_o;
        _snp_name_map=snp_name_map_o;
        LOGGER<<"-----------------------------------"<<std::endl;
    }
    
    for(c1=0; c1<chrs.size(); c1++){
        delete[] (grm_chrs[c1]);
        delete[] (geno_chrs[c1]);
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


void gcta::save_reml_state(std::string filename, bool no_adj_covar)
{
    int n = _n;
    int x_c = _X_c;
    int num_varcmp = static_cast<int>(_varcmp.size());
    int num_r_indx = static_cast<int>(_r_indx.size());

    std::ofstream outfile(filename, std::ios::binary);
    if(!outfile.is_open()) LOGGER.e(0, "cannot open the file ["+filename+"] to write.");

    // magic + dimensions
    outfile.write("REML", 4);
    outfile.write(reinterpret_cast<const char*>(&n),          sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&x_c),        sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&num_varcmp), sizeof(int));
    outfile.write(reinterpret_cast<const char*>(&num_r_indx), sizeof(int));

    // _Vi (row-major) — bulk write via a float buffer
    {
        size_t vi_floats = static_cast<size_t>(n) * n;
        std::vector<float> vi_buf(vi_floats);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                vi_buf[static_cast<size_t>(i) * n + j] = static_cast<float>(_Vi(i, j));
        outfile.write(reinterpret_cast<const char*>(vi_buf.data()),
                      static_cast<std::streamsize>(vi_floats * sizeof(float)));
    }

    // _b
    if(!no_adj_covar) {
        std::vector<float> b_buf(x_c);
        for(int i = 0; i < x_c; i++) b_buf[i] = static_cast<float>(_b(i));
        outfile.write(reinterpret_cast<const char*>(b_buf.data()), x_c * sizeof(float));
    }

    // _varcmp
    {
        std::vector<float> vc_buf(num_varcmp);
        for(int i = 0; i < num_varcmp; i++) vc_buf[i] = static_cast<float>(_varcmp[i]);
        outfile.write(reinterpret_cast<const char*>(vc_buf.data()), num_varcmp * sizeof(float));
    }

    // _r_indx
    outfile.write(reinterpret_cast<const char*>(_r_indx.data()), num_r_indx * sizeof(int));

    outfile.close();
    LOGGER << "Saved REML state (n=" << n << ", covariates=" << x_c
           << ", variance components=" << num_varcmp << ") to [" << filename << "]." << std::endl;
}

void gcta::load_reml_state(std::string filename, bool no_adj_covar)
{
    std::ifstream infile(filename, std::ios::binary);
    if(!infile.is_open()) LOGGER.e(0, "cannot open the file ["+filename+"] to read. Make sure you have run --mlma --save-reml first with matching --out prefix.");

    // magic
    char magic[4];
    infile.read(magic, 4);
    if(std::string_view(magic, 4) != "REML")
        LOGGER.e(0, "file ["+filename+"] is not a valid REML state file.");

    // dimensions
    int n = 0, x_c = 0, num_varcmp = 0, num_r_indx = 0;
    infile.read(reinterpret_cast<char*>(&n),          sizeof(int));
    infile.read(reinterpret_cast<char*>(&x_c),        sizeof(int));
    infile.read(reinterpret_cast<char*>(&num_varcmp), sizeof(int));
    infile.read(reinterpret_cast<char*>(&num_r_indx), sizeof(int));

    if(n != _n)
        LOGGER.e(0, "sample size mismatch: REML state has n=" + std::to_string(n) + " but current dataset has n=" + std::to_string(_n));
    if(x_c != _X_c)
        LOGGER.e(0, "number of covariates mismatch: REML state has " + std::to_string(x_c) + " covariates but current dataset has " + std::to_string(_X_c));

    // _Vi (row-major) — bulk read
    _Vi.resize(n, n);
    {
        size_t vi_floats = static_cast<size_t>(n) * n;
        std::vector<float> vi_buf(vi_floats);
        infile.read(reinterpret_cast<char*>(vi_buf.data()),
                    static_cast<std::streamsize>(vi_floats * sizeof(float)));
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                _Vi(i, j) = vi_buf[static_cast<size_t>(i) * n + j];
    }

    // _b
    if(!no_adj_covar) {
        _b.resize(x_c);
        std::vector<float> b_buf(x_c);
        infile.read(reinterpret_cast<char*>(b_buf.data()), x_c * sizeof(float));
        for(int i = 0; i < x_c; i++) _b(i) = b_buf[i];
    }

    // _varcmp
    _varcmp.resize(num_varcmp);
    {
        std::vector<float> vc_buf(num_varcmp);
        infile.read(reinterpret_cast<char*>(vc_buf.data()), num_varcmp * sizeof(float));
        for(int i = 0; i < num_varcmp; i++) _varcmp[i] = vc_buf[i];
    }

    // _r_indx
    _r_indx.resize(num_r_indx);
    infile.read(reinterpret_cast<char*>(_r_indx.data()), num_r_indx * sizeof(int));

    infile.close();
    LOGGER << "Loaded REML state from [" << filename << "]: n=" << n
           << ", covariates=" << x_c << ", variance components=" << num_varcmp << std::endl;
}
