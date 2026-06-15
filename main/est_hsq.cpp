/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for REML analysis
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"
#include "mem.hpp"
#include "StatFunc.h"
#include <random>

void gcta::set_reml_diag_mul(double value){
    _reml_diag_mul = value;
}

void gcta::set_reml_diagV_adj(int method){
    _reml_diagV_adj = method;
}

void gcta::set_reml_force_inv()
{
    _reml_force_inv = true;
}

void gcta::set_log_pval(bool log_pval)
{
    _log_pval = log_pval;
}

void gcta::set_reml_allow_constrain_run(){
    _reml_allow_constrain_run = true;
}

void gcta::set_reml_force_converge()
{
    _reml_force_converge = true;
}

void gcta::set_cv_blup(bool cv_blup){
    _cv_blup = cv_blup;
}

void gcta::set_reml_no_converge()
{
    _reml_no_converge = true;
}

void gcta::set_reml_fixed_var()
{
    _reml_fixed_var = true;
}

void gcta::set_reml_mtd(int reml_mtd)
{
    _reml_mtd = reml_mtd;
}

void gcta::set_reml_inv_method(int method){
    _reml_inv_mtd = method;
}

void gcta::set_reml_trace_approx(bool on, int nprobes){
    _reml_trace_approx = on;
    _reml_trace_approx_nprobes = nprobes;
}

void gcta::set_reml_trace_power_iter(int q){
    _reml_trace_power_iter = q;
}
    

void gcta::read_phen(std::string phen_file, std::vector<std::string> &phen_ID, std::vector< std::vector<std::string> > &phen_buf, int mphen, int mphen2) {
    // Read phenotype data
    std::ifstream in_phen(phen_file.c_str());
    if (!in_phen) LOGGER.e(0, "cannot open the file [" + phen_file + "] to read.");

    int i = 0;
    std::vector<std::string> fid, pid, vs_buf;
    std::string str_buf, fid_buf, pid_buf;
    phen_ID.clear();
    phen_buf.clear();
    LOGGER << "Reading phenotypes from [" + phen_file + "]." << std::endl;
    std::getline(in_phen, str_buf);
    int phen_num = StrFunc::split_string(str_buf, vs_buf) - 2;
    if (phen_num <= 0) LOGGER.e(0, "no phenotype data is found.");
    if (phen_num > 1) LOGGER << "There are " << phen_num << " traits specified in the file [" + phen_file + "]." << std::endl;
    if (mphen > phen_num) {
        std::stringstream errmsg;
        errmsg << "can not find the " << mphen << "th trait in the file [" + phen_file + "].";
        LOGGER.e(0, errmsg.str());
    }
    if (_bivar_reml && mphen2 > phen_num) {
        std::stringstream errmsg;
        errmsg << "can not find the " << mphen2 << "th trait in the file [" + phen_file + "].";
        LOGGER.e(0, errmsg.str());
    }
    if (_bivar_reml) LOGGER << "Traits " << mphen << " and " << mphen2 << " are included in the bivariate analysis." << std::endl;
    else {
        if (phen_num > 1) LOGGER << "Trait #" << mphen << " is included for analysis." << std::endl;
    }
    in_phen.seekg(std::ios::beg);
    mphen--;
    mphen2--;
    int line = 1;
    while (in_phen) {
        line++;
        in_phen >> fid_buf;
        if (in_phen.eof()) break;
        in_phen >> pid_buf;
        std::getline(in_phen, str_buf);
        if (StrFunc::split_string(str_buf, vs_buf) != phen_num) {
            std::stringstream errmsg;
            errmsg << vs_buf.size() - phen_num << " phenotype values are missing in line #" << line << " in the file [" + phen_file + "]";
            LOGGER.e(0, errmsg.str());
        }
        if (_bivar_reml) {
            if ((vs_buf[mphen] == "-9" || vs_buf[mphen] == "NA") && (vs_buf[mphen2] == "-9" || vs_buf[mphen2] == "NA")) continue;
        } else {
            if (vs_buf[mphen] == "-9" || vs_buf[mphen] == "NA") continue;
        }
        phen_ID.push_back(fid_buf + ":" + pid_buf);
        fid.push_back(fid_buf);
        pid.push_back(pid_buf);
        phen_buf.push_back(vs_buf);
    }
    in_phen.close();
    LOGGER << "Non-missing phenotypes of " << phen_buf.size() << " individuals are included from [" + phen_file + "]." << std::endl;

    if (_id_map.empty()) {
        _fid = fid;
        _pid = pid;
        _indi_num = _fid.size();
        init_keep();
    }
}

int gcta::read_covar(std::string covar_file, std::vector<std::string> &covar_ID, std::vector< std::vector<std::string> > &covar, bool qcovar_flag) {
    // Read covariate data
    std::ifstream in_covar(covar_file.c_str());
    if (!in_covar) LOGGER.e(0, "cannot open the file [" + covar_file + "] to read.");

    int i = 0, covar_num = 0;
    std::string str_buf, id_buf;
    std::vector<std::string> covar_buf, vs_buf;
    covar_ID.clear();
    covar.clear();
    if (qcovar_flag) LOGGER << "Reading quantitative covariate(s) from [" + covar_file + "]." << std::endl;
    else LOGGER << "Reading discrete covariate(s) from [" + covar_file + "]." << std::endl;
    covar_num = read_fac(in_covar, covar_ID, covar);
    if (qcovar_flag) LOGGER << covar_num << " quantitative covariate(s) of " << covar_ID.size() << " individuals are included from [" + covar_file + "]." << std::endl;
    else LOGGER << covar_num << " discrete covariate(s) of " << covar_ID.size() << " individuals are included from [" + covar_file + "]." << std::endl;

    return covar_num;
}

int gcta::read_fac(std::ifstream &ifstrm, std::vector<std::string> &ID, std::vector< std::vector<std::string> > &fac) {
    int i = 0, line = 0, fac_num = 0, prev_fac_num = 0;
    std::string str_buf, id_buf;
    std::vector<std::string> vs_buf;
    while (ifstrm) {
        ifstrm >> str_buf;
        if (ifstrm.eof()) break;
        id_buf = str_buf;
        ifstrm >> str_buf;
        id_buf += ":" + str_buf;
        std::getline(ifstrm, str_buf);
        fac_num = StrFunc::split_string(str_buf, vs_buf);
        if (line > 0 && fac_num != prev_fac_num) LOGGER.e(0, "each row should have the same number of columns.\n" + id_buf + "\t" + str_buf);
        line++;
        prev_fac_num = fac_num;
        bool continue_flag = false;
        for (i = 0; i < fac_num; i++) {
            if (vs_buf[i] == "-9" || vs_buf[i] == "NA") continue_flag = true;
        }
        if (continue_flag) continue;
        ID.push_back(id_buf);
        fac.push_back(vs_buf);
    }
    ifstrm.close();
    return fac_num;
}

int gcta::read_GE(std::string GE_file, std::vector<std::string> &GE_ID, std::vector< std::vector<std::string> > &GE, bool qGE_flag) {
    // Read phenotype data
    std::ifstream in_GE(GE_file.c_str());
    if (!in_GE) LOGGER.e(0, "cannot open the file [" + GE_file + "] to read.");

    std::string str_buf, id_buf;
    std::vector<std::string> vs_buf;
    GE_ID.clear();
    GE.clear();
    std::string env = "environmental";
    if (qGE_flag == true) env = "continuous " + env;
    else env = "categorical " + env;
    LOGGER << "Reading " << env << " factor(s) for the analysis of GE interaction from [" + GE_file + "]." << std::endl;
    int GE_num = read_fac(in_GE, GE_ID, GE);
    if (GE_num == 0) LOGGER.e(0, "no " + env + " factor is specified. Please check the format of the file: " + GE_file + ".");
    LOGGER << GE_num << " " << env << " factor(s) for " << GE_ID.size() << " individuals are included from [" + GE_file + "]." << std::endl;

    return GE_num;
}

void gcta::read_weight(std::string phen_file, std::vector<std::string> &phen_ID, std::vector<double> &weights) {
    std::ifstream in_phen(phen_file.c_str());
    if (!in_phen) LOGGER.e(0, "cannot open the weight [" + phen_file + "] to read.");

    phen_ID.clear();
    weights.clear();
    LOGGER << "Reading weights from [" + phen_file + "]." << std::endl;

    std::string str_buf, fid_buf, pid_buf;
    std::vector<std::string> vs_buf;

    // Parse header to determine column count
    std::getline(in_phen, str_buf);
    int phen_num = StrFunc::split_string(str_buf, vs_buf) - 2;
    if (phen_num <= 0) LOGGER.e(0, "no weight data is found.");
    if (phen_num > 1)
        LOGGER << "There are " << phen_num << " weights in [" + phen_file + "], only the first is used." << std::endl;

    // seekg(0, beg) — seekg(ios::beg) is seekg(0) by coincidence, not by spec
    in_phen.seekg(0, std::ios::beg);

    int line = 1;
    while (in_phen >> fid_buf >> pid_buf && std::getline(in_phen, str_buf)) {
        ++line;
        if (StrFunc::split_string(str_buf, vs_buf) != phen_num) {
            std::stringstream errmsg;
            errmsg << vs_buf.size() - phen_num << " weight values missing in line #" << line
                   << " in [" + phen_file + "]";
            LOGGER.e(0, errmsg.str());
        }

        const char* temp_str = vs_buf[0].c_str();
        char* pEnd;
        double temp_double = strtod(temp_str, &pEnd);
        if (std::isfinite(temp_double) && strlen(temp_str) == static_cast<size_t>(pEnd - temp_str)) {
            phen_ID.push_back(fid_buf + ":" + pid_buf);
            weights.push_back(temp_double);
        } else {
            LOGGER.w(0, "ignored line #" + std::to_string(line) + ".");
        }
    }

    LOGGER << "Non-missing weights of " << phen_ID.size()
           << " individuals included from [" + phen_file + "]." << std::endl;
}


void gcta::fit_reml(std::string grm_file, std::string phen_file, std::string qcovar_file, std::string covar_file, std::string qGE_file, std::string GE_file, std::string keep_indi_file, std::string remove_indi_file, std::string sex_file, int mphen, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool m_grm_flag, bool pred_rand_eff, bool est_fix_eff, bool est_fix_eff_var, int reml_mtd, int MaxIter, std::vector<double> reml_priors, std::vector<double> reml_priors_var, std::vector<int> drop, bool no_lrt, double prevalence, bool no_constrain, bool mlmassoc, bool within_family, bool reml_bending, bool reml_diag_one, std::string weight_file) {
    _within_family = within_family;
    _reml_mtd = reml_mtd;
    _reml_max_iter = MaxIter;
    _reml_diag_one = reml_diag_one;
    bool grm_flag = (!grm_file.empty());
    bool qcovar_flag = (!qcovar_file.empty());
    bool covar_flag = (!covar_file.empty());
    bool GE_flag = (!GE_file.empty());
    bool qGE_flag = (!qGE_file.empty());
    if (m_grm_flag) grm_flag = false;

    // Read data
    std::stringstream errmsg;
    int qcovar_num = 0, covar_num = 0, qE_fac_num = 0, E_fac_num = 0;
    std::vector<std::string> phen_ID, qcovar_ID, covar_ID, qGE_ID, GE_ID, grm_id, grm_files;
    std::vector< std::vector<std::string> > phen_buf, qcovar, covar, GE, qGE; // save individuals by column

    if (grm_flag) {
        read_grm(grm_file, grm_id, true, false, !(adj_grm_fac > -1.0));
        update_id_map_kp(grm_id, _id_map, _keep);
        grm_files.push_back(grm_file);
    } 
    else if (m_grm_flag) {
        read_grm_filenames(grm_file, grm_files, false);
        for (int i = 0; i < grm_files.size(); i++) {
            read_grm(grm_files[i], grm_id, false, true, !(adj_grm_fac > -1.0));
            update_id_map_kp(grm_id, _id_map, _keep);
        }
    }
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if (qcovar_flag) {
        qcovar_num = read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if (covar_flag) {
        covar_num = read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    if (qGE_flag) {
        qE_fac_num = read_GE(qGE_file, qGE_ID, qGE, true);
        update_id_map_kp(qGE_ID, _id_map, _keep);
    }
    if (GE_flag) {
        E_fac_num = read_GE(GE_file, GE_ID, GE, false);
        update_id_map_kp(GE_ID, _id_map, _keep);
    }

    std::vector<std::string> weight_ID;
    std::vector<double> weights;
    if(!weight_file.empty()){
        read_weight(weight_file, weight_ID, weights);
        update_id_map_kp(weight_ID, _id_map, _keep);
    }

    if (!mlmassoc) {
        if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
        if (!remove_indi_file.empty()) remove_indi(remove_indi_file);
    }
    if (grm_flag) {
        if (grm_cutoff>-1.0) rm_cor_indi(grm_cutoff);
        if (!sex_file.empty()) update_sex(sex_file);
        if (adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
        if (dosage_compen>-1) dc(dosage_compen);
        _grm_N.resize(0, 0);
    }

    std::vector<std::string> uni_id;
    std::map<std::string, int> uni_id_map;
    make_uni_id(uni_id, uni_id_map);
    _n = _keep.size();
    if (_n < 1) LOGGER.e(0, "no individual is in common among the input files.");

    // construct model terms
    _y.setZero(_n);
    for (int i = 0; i < phen_ID.size(); i++) {
        auto iter = uni_id_map.find(phen_ID[i]);
        if (iter == uni_id_map.end()) continue;
        _y[iter->second] = atof(phen_buf[i][mphen - 1].c_str());
    }

    // case-control
    _ncase = 0.0;
    _flag_CC = check_case_control(_ncase, _y);
    LOGGER << std::endl;
    if (_flag_CC) {
        //if(mlmassoc) LOGGER.e(0, "the option --mlm-assoc is valid for quantitative traits only.");
        if (!_bivar_reml) {
            if (!mlmassoc && prevalence<-1) LOGGER.i(0, "you can specify the disease prevalence by the option --prevalence so that GCTA can transform the variance explained to the underlying liability scale.", "Note:");
        } else {
            LOGGER.i(0, "you can specify the prevalences of the two diseases by the option --reml-bivar-prevalence so that GCTA can transform the estimates of variance explained to the underlying liability scale.", "Note:");
        }
    }


    int pos = 0;
    _r_indx.clear();
    std::vector<int> kp;
    if (grm_flag) {
        for (int i = 0; i < 1 + qE_fac_num + E_fac_num + 1; i++) _r_indx.push_back(i);
        if (!no_lrt) drop_comp(drop);
        _A.resize(_r_indx.size());
        if (mlmassoc) StrFunc::match(uni_id, grm_id, kp);
        else kp = _keep;
        (_A[0]) = eigenMatrix::Zero(_n, _n);

        {
            eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
            Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
            (_A[0]) = grm_sym(kp_idx, kp_idx);
        }
        if (_reml_diag_one) {
            double diag_mean = (_A[0]).diagonal().mean();
            LOGGER << "Mean of diagonal elements of the GRM = " << diag_mean << std::endl;
            //pragma omp parallel for private(j)
            {
                eigenVector diag_inv = (_A[0]).diagonal().cwiseInverse();
                (_A[0]) = ((_A[0]).array().colwise() / diag_inv.array()).matrix();
                (_A[0]).triangularView<Eigen::StrictlyUpper>() =
                    (_A[0]).triangularView<Eigen::StrictlyLower>().transpose();
            }
        }

        pos++;
        _grm.resize(0, 0);
    } else if (m_grm_flag) {
        if (!sex_file.empty()) update_sex(sex_file);
        for (int i = 0; i < (1 + qE_fac_num + E_fac_num) * grm_files.size() + 1; i++) _r_indx.push_back(i);
        if (!no_lrt) drop_comp(drop);
        _A.resize(_r_indx.size());
        std::string prev_file = grm_files[0];
        std::vector<std::string> prev_grm_id(grm_id);
        LOGGER << "There are " << grm_files.size() << " GRM file names specified in the file [" + grm_file + "]." << std::endl;
        for (int i = 0; i < grm_files.size(); i++, pos++) {
            LOGGER << "Reading the GRM from the " << i + 1 << "th file ..." << std::endl;
            read_grm(grm_files[i], grm_id, true, false, !(adj_grm_fac > -1.0));
            if (adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
            if (dosage_compen>-1) dc(dosage_compen);
            StrFunc::match(uni_id, grm_id, kp);
            (_A[pos]) = eigenMatrix::Zero(_n, _n);

            {
                eigenMatrix grm_sym(_grm.selfadjointView<Eigen::Lower>());
                Eigen::Map<const Eigen::VectorXi> kp_idx(kp.data(), _n);
                (_A[pos]) = grm_sym(kp_idx, kp_idx);
            }

            if (_reml_diag_one) {
                double diag_mean = (_A[pos]).diagonal().mean();
                LOGGER << "Mean of diagonal elements of the GRM = " << diag_mean << std::endl;
                {
                    eigenVector diag_inv = (_A[pos]).diagonal().cwiseInverse();
                    (_A[pos]) = ((_A[pos]).array().colwise() / diag_inv.array()).matrix();
                    (_A[pos]).triangularView<Eigen::StrictlyUpper>() =
                        (_A[pos]).triangularView<Eigen::StrictlyLower>().transpose();
                }
            }

            prev_file = grm_files[i];
            prev_grm_id = grm_id;
        }
        _grm_N.resize(0, 0);
        _grm.resize(0, 0);
    }
    else {
        _r_indx.push_back(0);
        _A.resize(_r_indx.size());
    }
    // Do NOT unconditionally allocate eigenMatrix::Identity here: it can be large allocation for a matrix whose only effect is σ²_e * I.  Convention:
    // _A[ci].size() == 0 means "identity component"; every hot path (calcu_Vi,
    // ai_reml, calcu_tr_PA_hutchpp, etc.) handles this case explicitly.  The full
    // n×n matrix is only allocated when a weight file overrides the diagonal.
    if(!weight_file.empty()){
        _A[_r_indx.size() - 1] = eigenMatrix::Identity(_n, _n);
        // contruct weight
        Eigen::VectorXd v_weight(_n);
        for (int i = 0; i < weight_ID.size(); i++) {
            auto iter = uni_id_map.find(weight_ID[i]);
            if (iter == uni_id_map.end()) continue;
            v_weight(iter->second) = weights[i];
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

    // GE interaction
    std::vector<eigenMatrix> E_float(E_fac_num);
    eigenMatrix qE_float, mbuf;
    if (qGE_flag) {
        qE_float.resize(_n, qE_fac_num);
        for (int i = 0; i < qGE_ID.size(); i++) {
            auto iter = uni_id_map.find(qGE_ID[i]);
            if (iter == uni_id_map.end()) continue;
            for (int j = 0; j < qE_fac_num; j++) qE_float(iter->second, j) = atof(qGE[i][j].c_str());
        }
        for (int j = 0; j < qE_fac_num; j++) {
            mbuf = ((qE_float.block(0, j, _n, 1))*(qE_float.block(0, j, _n, 1)).transpose());
            for (int i = 0; i < grm_files.size(); i++, pos++) (_A[pos]) = (_A[i]).array() * mbuf.array();
        }
    }
    if (GE_flag) {
        std::vector< std::vector<std::string> > E_str(E_fac_num);
        for (int i = 0; i < E_fac_num; i++) E_str[i].resize(_n);
        for (int i = 0; i < GE_ID.size(); i++) {
            auto iter = uni_id_map.find(GE_ID[i]);
            if (iter != uni_id_map.end()) {
                for (int j = 0; j < E_fac_num; j++) E_str[j][iter->second] = GE[i][j];
            }
        }
        for (int j = 0; j < E_fac_num; j++) {
            std::stringstream errmsg;
            errmsg << "too many classes for the " << j + 1 << "th environmental factor. \nPlease make sure you input a discrete variable as the environmental factor.";
            std::string errmsg1 = errmsg.str();
            errmsg.str("");
            errmsg << "the " << j + 1 << "th environmental factor has only one class.";
            std::string errmsg2 = errmsg.str();
            coeff_mat(E_str[j], E_float[j], errmsg1, errmsg2);
            mbuf = ((E_float[j])*(E_float[j]).transpose());
            for (int i = 0; i < grm_files.size(); i++, pos++) (_A[pos]) = (_A[i]).array() * mbuf.array();
        }
    }

    // construct X matrix
    construct_X(_n, uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);

    // names of variance component
    for (int i = 0; i < grm_files.size(); i++) {
        std::stringstream strstrm;
        if (grm_files.size() == 1) strstrm << "";
        else strstrm << i + 1;
        _var_name.push_back("V(G" + strstrm.str() + ")");
        _hsq_name.push_back("V(G" + strstrm.str() + ")/Vp");
    }
    for (int j = 0; j < qE_fac_num; j++) {
        for (int i = 0; i < grm_files.size(); i++) {
            std::stringstream strstrm1, strstrm2;
            if (grm_files.size() == 1) strstrm1 << "";
            else strstrm1 << i + 1;
            if (qE_fac_num == 1) strstrm2 << "";
            else strstrm2 << j + 1;
            _var_name.push_back("V(G" + strstrm1.str() + "xqE" + strstrm2.str() + ")");
            _hsq_name.push_back("V(G" + strstrm1.str() + "xqE" + strstrm2.str() + ")" + "/Vp");
        }
    }
    for (int j = 0; j < E_fac_num; j++) {
        for (int i = 0; i < grm_files.size(); i++) {
            std::stringstream strstrm1, strstrm2;
            if (grm_files.size() == 1) strstrm1 << "";
            else strstrm1 << i + 1;
            if (E_fac_num == 1) strstrm2 << "";
            else strstrm2 << j + 1;
            _var_name.push_back("V(G" + strstrm1.str() + "xE" + strstrm2.str() + ")");
            _hsq_name.push_back("V(G" + strstrm1.str() + "xE" + strstrm2.str() + ")" + "/Vp");
        }
    }
    _var_name.push_back("V(e)");

    LOGGER << _n << " individuals are in common in these files." << std::endl;

    // within family
    if (_within_family) detect_family();

    // bending
    if (reml_bending) bend_A();
    //LOGGER << "Prepare time: " << LOGGER.tp("main") << std::endl;

    // run REML algorithm
    if (_woodbury_rank > 0)
        compute_woodbury_basis(_woodbury_rank, 0.0, 0);
    else if (_woodbury_rank < 0)
        compute_woodbury_basis(0, _woodbury_buffer_factor, _woodbury_k_max);
    reml(pred_rand_eff, est_fix_eff, est_fix_eff_var, reml_priors, reml_priors_var, prevalence, -2.0, no_constrain, no_lrt, mlmassoc);
}

void gcta::drop_comp(std::vector<int> &drop) {
    int i = 0;
    std::stringstream errmsg;
    _r_indx_drop = _r_indx;
    std::stable_sort(drop.begin(), drop.end());
    drop.erase(unique(drop.begin(), drop.end()), drop.end());
    for (i = drop.size() - 1; i >= 0; i--) {
        if (drop[i] < 1 || drop[i] > _r_indx.size() - 1) {
            errmsg << "there " << (_r_indx.size() > 2 ? "are" : "is") << " only " << _r_indx.size() - 1 << " genetic variance component in the model. You can't drop the " << drop[i] << "-th component.";
            LOGGER.e(0, errmsg.str());
        }
        _r_indx_drop.erase(_r_indx_drop.begin() + drop[i] - 1);
    }
    if (_r_indx.size() == _r_indx_drop.size()) LOGGER.e(0, "no component has been dropped from the model. Please check the --reml-lrt option.");
}

void gcta::construct_X(int n, std::map<std::string, int> &uni_id_map, bool qcovar_flag, int qcovar_num, std::vector<std::string> &qcovar_ID, std::vector< std::vector<std::string> > &qcovar, bool covar_flag, int covar_num, std::vector<std::string> &covar_ID, std::vector< std::vector<std::string> > &covar, std::vector<eigenMatrix> &E_float, eigenMatrix &qE_float) {
    int i = 0, j = 0;
    std::stringstream errmsg;

    _X_c = 1;
    // quantitative covariates
    eigenMatrix X_q;
    if (qcovar_flag) {
        X_q.resize(n, qcovar_num);
        for (i = 0; i < qcovar_ID.size(); i++) {
            auto iter = uni_id_map.find(qcovar_ID[i]);
            if (iter == uni_id_map.end()) continue;
            for (j = 0; j < qcovar_num; j++) X_q(iter->second, j) = atof(qcovar[i][j].c_str());
        }
        if (qcovar_num == 0) LOGGER.e(0, "no quantitative covariate is found.");
        LOGGER << qcovar_num << " quantitative variable(s) included as covariate(s)." << std::endl;
        _X_c += qcovar_num;
    }
    // discrete covariates
    std::vector<eigenMatrix> X_d;
    if (covar_flag) {
        std::vector< std::vector<std::string> > covar_tmp(covar_num);
        for (i = 0; i < covar_num; i++) covar_tmp[i].resize(n);
        for (i = 0; i < covar_ID.size(); i++) {
            auto iter = uni_id_map.find(covar_ID[i]);
            if (iter == uni_id_map.end()) continue;
            for (j = 0; j < covar_num; j++) covar_tmp[j][iter->second] = covar[i][j];
        }
        LOGGER << covar_num << " discrete variable(s) included as covariate(s)." << std::endl;
        X_d.resize(covar_num);
        for (i = 0; i < covar_num; i++) {
            std::stringstream errmsg;
            errmsg << "too many classes for the " << i + 1 << "th discrete variable. \nPlease use the --qcovar if it is a quantitative covariate.";
            std::string errmsg1 = errmsg.str();
            errmsg.str("");
            errmsg << "the " << i + 1 << "th discrete variable has only one class.";
            std::string errmsg2 = errmsg.str();
            coeff_mat(covar_tmp[i], X_d[i], errmsg1, errmsg2);
            _X_c += (X_d[i]).cols() - 1;
        }
    }

    // E factor
    _X_c += qE_float.cols();
    for (i = 0; i < E_float.size(); i++) _X_c += (E_float[i]).cols() - 1;

    // Construct _X
    int col = 0;
    _X.resize(n, _X_c);
    _X.block(0, col, n, 1) = eigenMatrix::Ones(n, 1);
    col++;
    if (qcovar_flag) {
        _X.block(0, col, n, X_q.cols()) = X_q;
        col += X_q.cols();
    }
    for (i = 0; i < X_d.size(); i++) {
        _X.block(0, col, n, (X_d[i]).cols() - 1) = (X_d[i]).block(0, 1, n, (X_d[i]).cols() - 1);
        col += (X_d[i]).cols() - 1;
    }
    if (qE_float.cols() > 0) {
        _X.block(0, col, n, qE_float.cols()) = qE_float;
        col += qE_float.cols();
    }
    for (i = 0; i < E_float.size(); i++) {
        _X.block(0, col, n, (E_float[i]).cols() - 1) = (E_float[i]).block(0, 1, n, (E_float[i]).cols() - 1);
        col += (E_float[i]).cols() - 1;
    }
}

void gcta::coeff_mat(const std::vector<std::string> &vec, eigenMatrix &coeff_mat, std::string errmsg1, std::string errmsg2) {
    std::vector<std::string> value(vec);
    std::stable_sort(value.begin(), value.end());
    value.erase(unique(value.begin(), value.end()), value.end());
    if (value.size() > 0.5 * vec.size()) LOGGER.e(0, errmsg1); // LOGGER.e(0, "too many classes for the environmental factor. \nPlease make sure you input a discrete variable as the environmental factor.");
    if (value.size() == 1) LOGGER.e(0, errmsg2); //LOGGER.e(0, "the environmental factor should have at least two classes.");

    int i = 0, j = 0, row_num = vec.size(), column_num = value.size();
    std::map<std::string, int> val_map;
    for (i = 0; i < value.size(); i++) val_map.insert(std::pair<std::string, int>(value[i], i));

    coeff_mat.resize(row_num, column_num);
    coeff_mat.setZero(row_num, column_num);
    for (i = 0; i < row_num; i++) {
        auto iter = val_map.find(vec[i]);
        coeff_mat(i, iter->second) = 1.0;
    }
}

bool gcta::check_case_control(double &ncase, eigenVector &y) {
    int n = y.size();
    double case_num = 0.0;
    std::vector<double> value(n);
    for (int i = 0; i < n; i++) value[i] = y(i);
    std::stable_sort(value.begin(), value.end());
    value.erase(unique(value.begin(), value.end()), value.end());
    if (value.size() == 2) {
        if (CommFunc::FloatEqual(value[0], 0.0) && CommFunc::FloatEqual(value[1], 1.0)) case_num = y.sum();
        else if (CommFunc::FloatEqual(value[0], 1.0) && CommFunc::FloatEqual(value[1], 2.0)) case_num = (y.sum() - n);
        if (!_bivar_reml) LOGGER  << "Assuming a disease phenotype for a case-control study: ";
        LOGGER << (int) case_num << " cases and " << (int) (n - case_num) << " controls ";
        ncase = case_num / (double) n;
        return true;
    } else if (value.size() < 2) {
        LOGGER << "Only " << value.size() << " type(s) of phenotype values observed: ";
        for(int i = 0; i < value.size(); i++){
            LOGGER << value[i] << "\t";
        }
        LOGGER << std::endl;
        LOGGER.e(0, "invalid phenotype. Please check the phenotype file.");
    }
    return false;
}

double gcta::transform_hsq_L(double P, double K, double hsq) {
    double t = StatFunc::qnorm(1.0 - K);
    double z = StatFunc::dnorm(t);
    double C = (K * (1 - K) / (z * z))*(K * (1 - K) / (P * (1 - P)));
    return (hsq * C);
}

int gcta::constrain_varcmp(eigenVector &varcmp) {
    int pos = 0;
    double delta = 0.0, constr_scale = 1e-6;
    int i = 0, num = 0;
    std::vector<int> constrain(_r_indx.size());

    if (_bivar_reml) {
        for (i = 0, num = 0; i < _bivar_pos[0].size(); i++) {
            pos = _bivar_pos[0][i];
            if (varcmp[pos] < 0) {
                delta += _y_Ssq * constr_scale - varcmp[pos];
                varcmp[pos] = _y_Ssq * constr_scale;
                constrain[i] = 1;
                num++;
            }
        }
        delta /= (_bivar_pos[0].size() - num);
        for (i = 0; i < _bivar_pos[0].size(); i++) {
            pos = _bivar_pos[0][i];
            if (constrain[pos] < 1 && varcmp[pos] > delta) varcmp[pos] -= delta;
        }

        for (i = 0, num = 0; i < _bivar_pos[1].size(); i++) {
            pos = _bivar_pos[1][i];
            if (varcmp[pos] < 0) {
                delta += _y_Ssq * constr_scale - varcmp[pos];
                varcmp[pos] = _y_Ssq * constr_scale;
                constrain[i] = 1;
                num++;
            }
        }
        delta /= (_bivar_pos[1].size() - num);
        for (i = 0; i < _bivar_pos[1].size(); i++) {
            pos = _bivar_pos[1][i];
            if (constrain[pos] < 1 && varcmp[pos] > delta) varcmp[pos] -= delta;
        }

        for (i = 0, num = 0; i < constrain.size(); i++) {
            if (constrain[i] == 1) num++;

        }
        return num;
    }

    for (i = 0; i < _r_indx.size(); i++) {
        if (varcmp[i] < 0) {
            delta += _y_Ssq * constr_scale - varcmp[i];
            varcmp[i] = _y_Ssq * constr_scale;
            constrain[i] = 1;
            num++;
        }
    }
    delta /= (_r_indx.size() - num);
    for (i = 0; i < _r_indx.size(); i++) {
        if (constrain[i] < 1 && varcmp[i] > delta) varcmp[i] -= delta;
    }

    return num;
}

void gcta::reml(bool pred_rand_eff, bool est_fix_eff, bool est_fix_eff_var, std::vector<double> &reml_priors, std::vector<double> &reml_priors_var, double prevalence, double prevalence2, bool no_constrain, bool no_lrt, bool mlmassoc)
{
    int i = 0, j = 0, k = 0;

    // Initialize variance component
    // 0: AI; 1: Fisher-scoring; 2: EM
    std::stringstream errmsg;
    double d_buf = 0.0;
    eigenVector y_tmp = _y.array() - _y.mean();
    if (!_bivar_reml) {
        _y_Ssq = y_tmp.squaredNorm() / (_n - 1.0);
        if (!(fabs(_y_Ssq) < 1e30)) LOGGER.e(0, "the phenotypic variance is infinite. Please check the missing data in your phenotype file. Missing values should be represented by \"NA\" or \"-9\".");
    }
    bool reml_priors_flag = !reml_priors.empty(), reml_priors_var_flag = !reml_priors_var.empty();
    if (reml_priors_flag && reml_priors.size() < _r_indx.size() - 1) {
        errmsg << "in option --reml-priors. There are " << _r_indx.size() << " variance components. At least " << _r_indx.size() - 1 << " prior values should be specified.";
        LOGGER.e(0, errmsg.str());
    }
    if (reml_priors_var_flag && reml_priors_var.size() < _r_indx.size() - 1) {
        errmsg << "in option --reml-priors-var. There are " << _r_indx.size() << " variance components. At least " << _r_indx.size() - 1 << " prior values should be specified.";
        LOGGER.e(0, errmsg.str());
    }

    LOGGER << "\nPerforming " << (_bivar_reml ? "bivariate" : "") << " REML analysis ... (Note: may take hours depending on sample size)." << std::endl;
    if (_n < 10) LOGGER.e(0, "sample size is too small.");
    LOGGER << _n << " observations, " << _X_c << " fixed effect(s), and " << _r_indx.size() << " variance component(s)(including residual variance)." << std::endl;
    eigenMatrix Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(_r_indx.size(), _r_indx.size());
    eigenVector Py(_n), varcmp;
    init_varcomp(reml_priors_var, reml_priors, varcmp);
    //LOGGER << "REML begin reml iteration" << std::endl;
    double lgL = reml_iteration(Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, reml_priors_var_flag | reml_priors_flag, no_constrain);
    eigenMatrix u;
    if (pred_rand_eff) {
        u.resize(_n, _r_indx.size());
        for (i = 0; i < _r_indx.size(); i++) {
            if (_bivar_reml || _within_family)(u.col(i)) = (((_Asp[_r_indx[i]]) * Py) * varcmp[i]);
            else if (_A[_r_indx[i]].size() == 0) u.col(i) = Py * varcmp[i];
            else (u.col(i)) = (((_A[_r_indx[i]]) * Py) * varcmp[i]);
        }
    }
    if (est_fix_eff) _b = Xt_Vi_X_i * (Vi_X.transpose() * _y);
    
    // calculate Hsq and SE
    double Vp = 0.0, Vp2 = 0.0, VarVp = 0.0, VarVp2 = 0.0, Vp_f = 0.0, VarVp_f = 0.0;
    std::vector<double> Hsq(_r_indx.size() - 1), VarHsq(_r_indx.size() - 1);
    calcu_Vp(Vp, Vp2, VarVp, VarVp2, varcmp, Hi);
    for (i = 0; i < Hsq.size(); i++) calcu_hsq(i, Vp, Vp2, VarVp, VarVp2, Hsq[i], VarHsq[i], varcmp, Hi);

    // calculate the logL for a reduce model
    double lgL_rdu_mdl = 0.0, LRT = 0.0;
    if (!no_lrt) {
        lgL_rdu_mdl = lgL_reduce_mdl(no_constrain);
        if (lgL_rdu_mdl > lgL + 1e-6) {
            LOGGER.w(0, "reduced-model logL is greater than full-model logL; possible numerical instability.");
        }
        LRT = 2.0 * (lgL - lgL_rdu_mdl);
        if (LRT < 0.0) LRT = 0.0;
    }

    // calcuate the logL given a rG in a bivariate analysis
    double lgL_fixed_rg = 0.0;
    if (_bivar_reml && !_fixed_rg_val.empty()) {
        lgL_fixed_rg = lgL_fix_rg(varcmp, no_constrain);
        LRT = 2.0 * (lgL - lgL_fixed_rg);
        if (LRT < 0.0) LRT = 0.0;
    }

    // In the non-approx mlmassoc path the REML loop releases _Vi after every
    // iteration (including the last) to save n² memory between iterations.
    // Use factorize_only=true so both exact and Hutch++ paths enter mlma_calcu_stat
    // with L in _Vi_L and _Vi_use_llt=true, allowing mlma_calcu_stat to skip the
    // dpotri call entirely (see comment there).
    // The approx path (_Vi_use_llt=true) stores L in _Vi_L and is unaffected.
    if (mlmassoc && !_reml_trace_approx && !_Vi_use_woodbury && !_bivar_reml && !_within_family) {
        double logdet_dummy = 0.0;
        int iter_dummy = 0;
        calcu_Vi(_Vi, varcmp, logdet_dummy, iter_dummy, /*factorize_only=*/true);
    }

if (mlmassoc) {
    _varcmp = eigenVector2Vector(varcmp);

    std::string reml_rst_file = _out + ".hsq";
    std::ofstream o_reml(reml_rst_file.c_str());
    if (!o_reml) LOGGER.e(0, "cannot open the file [" + reml_rst_file + "] to write.");
    o_reml << "Source\tVariance\tSE\n";
    for (int i = 0; i < _r_indx.size(); i++)
        o_reml << _var_name[i] << "\t" << varcmp[i] << "\t" << std::sqrt(Hi(i,i)) << "\n";
    o_reml << "Vp\t" << Vp << "\t" << std::sqrt(VarVp) << "\n";
    for (int i = 0; i < Hsq.size(); i++)
        o_reml << _hsq_name[i] << "\t" << Hsq[i] << "\t" << std::sqrt(VarHsq[i]) << "\n";
    o_reml << "n\t" << _n << "\n";

    return;
}

    // output results
    double sum_hsq = 0.0, var_sum_hsq = 0.0;
    if (!_bivar_reml && _r_indx.size() > 2) calcu_sum_hsq(Vp, VarVp, sum_hsq, var_sum_hsq, varcmp, Hi);
    LOGGER << "\nSummary result of REML analysis:" << std::endl;
    LOGGER << "Source\tVariance\tSE" << std::fixed << LOGGER.setprecision(6) << std::endl;
    for (i = 0; i < _r_indx.size(); i++) LOGGER << _var_name[i] << "\t" << varcmp[i] << "\t" << sqrt(Hi(i, i)) << std::endl;
    if (_bivar_reml) {
        LOGGER << "Vp_tr1\t" << Vp << "\t" << sqrt(VarVp) << std::endl;
        LOGGER << "Vp_tr2\t" << Vp2 << "\t" << sqrt(VarVp2) << std::endl;
        for (i = 0, j = 0; i < _bivar_pos[0].size() - 1; i++, j += 2) {
            LOGGER << _hsq_name[j] << "\t" << Hsq[_bivar_pos[0][i]] << "\t" << sqrt(VarHsq[_bivar_pos[0][i]]) << std::endl;
            LOGGER << _hsq_name[j + 1] << "\t" << Hsq[_bivar_pos[1][i]] << "\t" << sqrt(VarHsq[_bivar_pos[1][i]]) << std::endl;
        }
    } else {
        LOGGER << "Vp\t" << Vp << "\t" << sqrt(VarVp) << std::endl;
        for (i = 0; i < Hsq.size(); i++) LOGGER << _hsq_name[i] << "\t" << Hsq[i] << "\t" << sqrt(VarHsq[i]) << std::endl;
        if (_r_indx.size() > 2) LOGGER << "\nSum of V(G)/Vp\t" << sum_hsq << "\t" << sqrt(var_sum_hsq) << std::endl;
    }
    if ((_flag_CC && prevalence>-1) || (_flag_CC2 && prevalence2>-1)) {
        LOGGER << "The estimate of variance explained on the observed scale is transformed to that on the underlying liability scale:" << std::endl;
        if (_bivar_reml) {
            if (_flag_CC) LOGGER << "Proportion of cases in the sample = " << _ncase << " for trait #1; User-specified disease prevalence = " << prevalence << " for trait #1" << std::endl;
            if (_flag_CC2) LOGGER << "Proportion of cases in the sample = " << _ncase2 << " for trait #2; User-specified disease prevalence = " << prevalence2 << " for trait #2" << std::endl;
            for (i = 0, j = 0; i < _bivar_pos[0].size() - 1; i++, j += 2) {
                if (_flag_CC) LOGGER << _hsq_name[j] << "_L\t" << transform_hsq_L(_ncase, prevalence, Hsq[_bivar_pos[0][i]]) << "\t" << transform_hsq_L(_ncase, prevalence, sqrt(VarHsq[_bivar_pos[0][i]])) << std::endl;
                if (_flag_CC2) LOGGER << _hsq_name[j + 1] << "_L\t" << transform_hsq_L(_ncase2, prevalence2, Hsq[_bivar_pos[1][i]]) << "\t" << transform_hsq_L(_ncase2, prevalence2, sqrt(VarHsq[_bivar_pos[1][i]])) << std::endl;
            }
        } else {
            LOGGER << "(Proportion of cases in the sample = " << _ncase << "; User-specified disease prevalence = " << prevalence << ")" << std::endl;
            for (i = 0; i < Hsq.size(); i++) LOGGER << _hsq_name[i] << "_L\t" << transform_hsq_L(_ncase, prevalence, Hsq[i]) << "\t" << transform_hsq_L(_ncase, prevalence, sqrt(VarHsq[i])) << std::endl;
            if (_r_indx.size() > 2)  LOGGER << "\nSum of V(G)_L/Vp\t" << transform_hsq_L(_ncase, prevalence, sum_hsq) << "\t" << transform_hsq_L(_ncase, prevalence, sqrt(var_sum_hsq)) << std::endl;
        }
    }
    // output genetic correlation
    eigenVector rg, rg_var;
    std::vector<std::string> rg_name;
    if (_bivar_reml) {
        calcu_rg(varcmp, Hi, rg, rg_var, rg_name);
        for (i = 0; i < rg_name.size(); i++) {
            LOGGER << rg_name[i] << "\t" << rg[i] << "\t" << sqrt(rg_var[i]) << std::endl;
        }
    }
    if(!_reml_force_converge || !_reml_AI_not_invertible){
        LOGGER << "\nSampling variance/covariance of the estimates of variance components:" << std::endl;
        for (i = 0; i < _r_indx.size(); i++) {
            //for (j = 0; j < _r_indx.size(); j++) LOGGER << setiosflags(std::ios::scientific) << Hi(i, j) << "\t";
            for (j = 0; j < _r_indx.size(); j++) LOGGER << std::scientific << Hi(i, j) << "\t";
            LOGGER << std::endl;
        }
    }
    if (est_fix_eff) {
        LOGGER << std::endl;
        LOGGER << "Estimate" << (_X_c > 1 ? "s" : "") << " of fixed effect" << (_X_c > 1 ? "s" : "") << ":" << std::endl;
        LOGGER << "\nSource\tEstimate\tSE" << std::endl;
        for (i = 0; i < _X_c; i++) {
            if (i == 0) LOGGER << "mean\t";
            else LOGGER << "X_" << i + 1 << "\t";
            LOGGER << std::fixed << _b[i] << "\t" << sqrt(Xt_Vi_X_i(i, i)) << std::endl;
        }
    }

    if (est_fix_eff_var) {
       LOGGER << std::endl;
       LOGGER << "Variance-covariance matrix of the estimated fixed effect"<< (_X_c > 1 ? "s" : "") << ":" << std::endl;
        for (i = 0; i < _X_c; i++) {
            for (j = 0; j < _X_c; j++) {
               LOGGER << std::fixed << Xt_Vi_X_i(i, j) << "\t";
            }
            LOGGER << std::endl;
        }
    }

    // save summary result into a file
    std::string reml_rst_file = _out + ".hsq";
    std::ofstream o_reml(reml_rst_file.c_str());
    if (!o_reml) LOGGER.e(0, "cannot open the file [" + reml_rst_file + "] to write.");
    o_reml << "Source\tVariance\tSE" << std::setiosflags(std::ios::fixed) << std::setprecision(6) << std::endl;
    for (i = 0; i < _r_indx.size(); i++) o_reml << _var_name[i] << "\t" << varcmp[i] << "\t" << sqrt(Hi(i, i)) << std::endl;
    if (_bivar_reml) {
        o_reml << "Vp_tr1\t" << Vp << "\t" << sqrt(VarVp) << std::endl;
        o_reml << "Vp_tr2\t" << Vp2 << "\t" << sqrt(VarVp2) << std::endl;
        for (i = 0, j = 0; i < _bivar_pos[0].size() - 1; i++, j += 2) {
            o_reml << _hsq_name[j] << "\t" << Hsq[_bivar_pos[0][i]] << "\t" << sqrt(VarHsq[_bivar_pos[0][i]]) << std::endl;
            o_reml << _hsq_name[j + 1] << "\t" << Hsq[_bivar_pos[1][i]] << "\t" << sqrt(VarHsq[_bivar_pos[1][i]]) << std::endl;
        }
    } else {
        o_reml << "Vp\t" << Vp << "\t" << sqrt(VarVp) << std::endl;
        for (i = 0; i < Hsq.size(); i++) o_reml << _hsq_name[i] << "\t" << Hsq[i] << "\t" << sqrt(VarHsq[i]) << std::endl;
        if (_r_indx.size() > 2) o_reml << "\nSum of V(G)/Vp\t" << sum_hsq << "\t" << sqrt(var_sum_hsq) << std::endl;
    }
    if (_flag_CC && prevalence>-1) {
        o_reml << "The estimate of variance explained on the observed scale is transformed to that on the underlying scale:" << std::endl;
        if (_bivar_reml) {
            o_reml << "(Proportion of cases in the sample = " << _ncase << "; User-specified disease prevalence = " << prevalence << " for disease 1 and = " << prevalence2 << " for disease 2)" << std::endl;
            for (i = 0, j = 0; i < _bivar_pos[0].size() - 1; i++, j += 2) {
                o_reml << _hsq_name[j] << "_L\t" << transform_hsq_L(_ncase, prevalence, Hsq[_bivar_pos[0][i]]) << "\t" << transform_hsq_L(_ncase, prevalence, sqrt(VarHsq[_bivar_pos[0][i]])) << std::endl;
                o_reml << _hsq_name[j + 1] << "_L\t" << transform_hsq_L(_ncase2, prevalence2, Hsq[_bivar_pos[1][i]]) << "\t" << transform_hsq_L(_ncase2, prevalence2, sqrt(VarHsq[_bivar_pos[1][i]])) << std::endl;
            }
        } else {
            o_reml << "(Proportion of cases in the sample = " << _ncase << "; User-specified disease prevalence = " << prevalence << ")" << std::endl;
            for (i = 0; i < Hsq.size(); i++) o_reml << _hsq_name[i] << "_L\t" << transform_hsq_L(_ncase, prevalence, Hsq[i]) << "\t" << transform_hsq_L(_ncase, prevalence, sqrt(VarHsq[i])) << std::endl;
            if (_r_indx.size() > 2)  o_reml << "\nSum of V(G)_L/Vp\t" << transform_hsq_L(_ncase, prevalence, sum_hsq) << "\t" << transform_hsq_L(_ncase, prevalence, sqrt(var_sum_hsq)) << std::endl;
        }
    }
    if (_bivar_reml) {
        for (i = 0; i < rg_name.size(); i++) o_reml << rg_name[i] << "\t" << rg[i] << "\t" << sqrt(rg_var[i]) << std::endl;
    }
    o_reml << "logL\t" << std::setprecision(3) << lgL << std::endl;
    if (!no_lrt && _r_indx.size() - 1 > 0) {
        o_reml << "logL0\t" << std::setprecision(3) << lgL_rdu_mdl << std::endl;
        o_reml << "LRT\t" << std::setprecision(3) << LRT << std::endl;
        o_reml << "df\t" << std::setprecision(1) << _r_indx.size() - _r_indx_drop.size() << std::endl;
        //o_reml << "Pval\t" << std::setprecision(4) << setiosflags(std::ios::scientific) << 0.5 * StatFunc::chi_prob(_r_indx.size() - _r_indx_drop.size(), LRT) << setiosflags(std::ios::fixed) << std::endl;
        if (_log_pval) {
            o_reml << "log_Pval\t" << std::setprecision(6) << (std::log(0.5) + StatFunc::chi_prob(_r_indx.size() - _r_indx_drop.size(), LRT, true)) << std::endl;
        } else {
            o_reml << "Pval\t" << std::scientific << std::setprecision(4) << 0.5 * StatFunc::chi_prob(_r_indx.size() - _r_indx_drop.size(), LRT) << std::fixed << std::endl;
        }
    }
    if (_bivar_reml && !_fixed_rg_val.empty()) {
        o_reml << "logL0\t" << std::setprecision(3) << lgL_fixed_rg << " (when rG fixed at ";
        for (i = 0; i < _fixed_rg_val.size() - 1; i++) o_reml << _fixed_rg_val[i] << "\t";
        o_reml << _fixed_rg_val[_fixed_rg_val.size() - 1] << ")" << std::endl;
        o_reml << "LRT\t" << std::setprecision(3) << LRT << std::endl;
        o_reml << "df\t" << std::setprecision(1) << _fixed_rg_val.size() << std::endl;
        //o_reml << "Pval\t" << std::setprecision(4) << setiosflags(std::ios::scientific) << 0.5 * StatFunc::chi_prob(_fixed_rg_val.size(), LRT) << setiosflags(std::ios::fixed) << " (one-tailed test)" << std::endl;
        if (_log_pval) {
            o_reml << "log_Pval\t" << std::setprecision(6) << (std::log(0.5) + StatFunc::chi_prob(_fixed_rg_val.size(), LRT, true)) << " (one-tailed test)" << std::endl;
        } else {
            o_reml << "Pval\t" << std::scientific << std::setprecision(4) << 0.5 * StatFunc::chi_prob(_fixed_rg_val.size(), LRT) << std::fixed << " (one-tailed test)" << std::endl;
        }
    }
    o_reml << "n\t" << _n << std::endl;
    if (est_fix_eff) {
        o_reml << "\nFix_eff\tSE" << std::endl;
        for (i = 0; i < _X_c; i++) o_reml << std::setprecision(6) << _b[i] << "\t" << sqrt(Xt_Vi_X_i(i, i)) << std::endl;
        o_reml.close();
    }
    LOGGER << "\nSummary result of REML analysis has been saved in the file [" + reml_rst_file + "]." << std::endl;

    // save random effect to a file
    if (pred_rand_eff) {
        std::string rand_eff_file = _out + ".indi.blp";
        std::ofstream o_rand_eff(rand_eff_file.c_str());
        for (i = 0; i < _keep.size(); i++) {
            o_rand_eff << _fid[_keep[i]] << "\t" << _pid[_keep[i]] << "\t";
            for (j = 0; j < _r_indx.size(); j++) o_rand_eff << std::setprecision(6) << Py[i] * varcmp[j] << "\t" << u(i, j) << "\t";
            o_rand_eff << std::endl;
        }
        LOGGER << "\nBLUP solutions of the genetic effects for " << _keep.size() << " individuals have been saved in the file [" + rand_eff_file + "]." << std::endl;
    }

    if(_cv_blup){
        std::string rand_eff_file = _out + ".indi.cvblp";
        std::ofstream o_rand_eff(rand_eff_file.c_str());
        double logdet; 
        int iter;
        // Other function changed the _Vi and some variables, thus need recalculate
        calcu_Vi(_Vi, varcmp, logdet, iter);
        calcu_P(_Vi, Vi_X, Xt_Vi_X_i, _P);

        eigenMatrix xt_vi_x_i_xtvi = Xt_Vi_X_i * Vi_X.transpose();
        // same above
        //eigenMatrix xt_vi_x_i_xtvi = Xt_Vi_X_i * (_X.transpose() * _Vi);
        _b = xt_vi_x_i_xtvi * _y; 
        eigenMatrix Xproj = _X * xt_vi_x_i_xtvi;
        eigenVector diag_Xproj = Xproj.diagonal();
        eigenVector y_fix_eff = Xproj * _y;

        // For test BLUP
        //eigenVector yFeResid = _y - y_fix_eff;
        //eigenMatrix H = varcmp[0] * _A[_r_indx[0]] * _Vi;
        //eigenVector gBlup = H * yFeResid;

        eigenVector y_fix_eff_loo = (y_fix_eff - diag_Xproj.cwiseProduct(_y)).array() / (1.0 - diag_Xproj.array()).array();
        eigenVector y_tilde = _y - y_fix_eff_loo;
        eigenVector y_tilde_centered = y_tilde.array() - y_tilde.mean();
        // _Vi stores only the lower triangle; selfadjointView routes to dsymv.
        eigenVector bBlup_base(_n);
        bBlup_base.noalias() = _Vi.selfadjointView<Eigen::Lower>() * y_tilde_centered;

        int col_num = _r_indx.size();
        eigenMatrix bBlups(bBlup_base.size(),  col_num);
        eigenMatrix cvBlups(bBlup_base.size(), col_num);

        // use the raw formula to calculate cvBLUP
        //eigenVector gBlupFeLoo = H * y_tilde_centered;
        //eigenVector cvBlups1 = (gBlupFeLoo.array() - H.diagonal().array() * y_tilde_centered.array()) / (1.0 - H.diagonal().array());
        

        for(int i = 0; i < col_num; i++){
            bBlups.col(i) = bBlup_base * varcmp[i];
            int _y_size = _y.size();
            eigenVector gBlup_fix_eff_loo(_y_size);
            eigenVector diag_H(_y_size);
            if (_A[_r_indx[i]].size() == 0) {
                // Identity component: A*v = v; diag(A * V^{-1}) = diag(V^{-1}).
                gBlup_fix_eff_loo = bBlups.col(i);
                for (int j = 0; j < _y_size; j++) diag_H[j] = _Vi(j, j) * varcmp[i];
            } else {
                gBlup_fix_eff_loo.noalias() = _A[_r_indx[i]] * bBlups.col(i);
                for(int j = 0; j < _y_size; j++){
                    // _Vi stores only the lower triangle.  Vi.col(j) decomposes as:
                    //   rows k > j : lower triangle, contiguous in column-major
                    //   rows k < j : upper triangle -> use Vi(j,k) by symmetry (row access)
                    double s = _Vi(j, j) * _A[_r_indx[i]](j, j);
                    if (j > 0)
                        s += _Vi.row(j).head(j).dot(_A[_r_indx[i]].row(j).head(j));
                    if (j + 1 < _y_size)
                        s += _Vi.col(j).tail(_y_size - j - 1).dot(_A[_r_indx[i]].col(j).tail(_y_size - j - 1));
                    diag_H[j] = s * varcmp[i];
                }
            }
            cvBlups.col(i) = (gBlup_fix_eff_loo.array() - diag_H.array() * y_tilde_centered.array()).array() / (1.0 - diag_H.array()).array();
        }

        // output the effects
        for(int i = 0; i < _keep.size(); i++){
            o_rand_eff << _fid[_keep[i]] << "\t" << _pid[_keep[i]] << "\t";
            for(int j = 0; j < col_num; j++){
                o_rand_eff << std::setprecision(6) << bBlups(i, j) << "\t" << cvBlups(i, j) << "\t";
            }
            //o_rand_eff << gBlup(i) << "\t" << cvBlups1(i);

            o_rand_eff << std::endl;
        }
        LOGGER << "\ncvBLUP solutions of the genetic effects for " << _keep.size() << " individuals have been saved in the file [" + rand_eff_file + "]." << std::endl;
    }

}

void gcta::init_varcomp(std::vector<double> &reml_priors_var, std::vector<double> &reml_priors, eigenVector &varcmp) {
    int i = 0, pos = 0;
    double d_buf = 0.0;

    varcmp = eigenVector::Zero(_r_indx.size());
    if (_bivar_reml) {
        if (!reml_priors_var.empty()) {
            for (i = 0; i < _r_indx.size(); i++) varcmp[i] = reml_priors_var[i];
        } 
        else if (!reml_priors.empty()) {
            for (i = 0, d_buf = 0; i < _bivar_pos[0].size() - 1; i++) {
                pos = _bivar_pos[0][i];
                varcmp[pos] = reml_priors[pos] * _y_Ssq;
                d_buf += reml_priors[pos];
            }
            if (d_buf > 1.0) LOGGER.e(0, "\n  --reml-priors. The sum of all prior values for trait 1 should not exceed 1.0.");
            varcmp[_bivar_pos[0][_bivar_pos[0].size() - 1]] = (1.0 - d_buf) * _y_Ssq;
            for (i = 0, d_buf = 0; i < _bivar_pos[1].size() - 1; i++) {
                pos = _bivar_pos[1][i];
                varcmp[pos] = reml_priors[pos] * _y_Ssq;
                d_buf += reml_priors[pos];
            }
            if (d_buf > 1.0) LOGGER.e(0, "\n  --reml-priors. The sum of all prior values for trait 2 should not exceed 1.0.");
            varcmp[_bivar_pos[1][_bivar_pos[1].size() - 1]] = (1.0 - d_buf) * _y2_Ssq;
            for (i = 0; i < _bivar_pos[2].size(); i++) varcmp[_bivar_pos[2][i]] = reml_priors[_bivar_pos[2][i]] * sqrt(_y_Ssq * _y2_Ssq);
        }
        else {
            for (i = 0; i < _bivar_pos[0].size(); i++) varcmp[_bivar_pos[0][i]] = _y_Ssq / _bivar_pos[0].size();
            for (i = 0; i < _bivar_pos[1].size(); i++) varcmp[_bivar_pos[1][i]] = _y2_Ssq / _bivar_pos[1].size();
            for (i = 0; i < _bivar_pos[2].size(); i++) varcmp[_bivar_pos[2][i]] = 0.5 * sqrt(varcmp[_bivar_pos[0][i]] * varcmp[_bivar_pos[1][i]]);
        }

        return;
    }

    if (!reml_priors_var.empty()) {
        for (i = 0; i < _r_indx.size() - 1; i++) varcmp[i] = reml_priors_var[i];
        if (reml_priors_var.size() < _r_indx.size()) varcmp[_r_indx.size() - 1] = _y_Ssq - varcmp.sum();
        else varcmp[_r_indx.size() - 1] = reml_priors_var[_r_indx.size() - 1];
    }
    else if (!reml_priors.empty()) {
        for (i = 0, d_buf = 0; i < _r_indx.size() - 1; i++) {
            varcmp[i] = reml_priors[i] * _y_Ssq;
            d_buf += reml_priors[i];
        }
        if (d_buf > 1.0) LOGGER.e(0, "\n  --reml-priors. The sum of all prior values should not exceed 1.0.");
        varcmp[_r_indx.size() - 1] = (1.0 - d_buf) * _y_Ssq;
    }
    else varcmp.setConstant(_y_Ssq / (_r_indx.size()));

    // Woodbury path: replace the flat default with an HE regression estimate.
    // The flat init (V(G) = V(e) = var(y)/2) is far from typical solutions (h² ≈ 0.1–0.3),
    // causing the single mandatory EM step at iter=0 to land poorly and requiring many
    // more AI-REML iterations to converge.
    //
    // HE estimator (assumes mean-centred y, negligible bias for typical GCTA phenotypes):
    //   sg = (n·y'Ky − tr(K)·y'y) / (n·tr(K²) − tr(K)²)
    // All quantities are O(nk) from woodbury members already computed.
    //
    // Only applied when no explicit priors were provided and the model is a standard
    // single-GRM univariate REML (not bivariate, not within-family).
    if (_Vi_use_woodbury && !_bivar_reml && !_within_family
            && (int)_r_indx.size() == 2
            && reml_priors.empty() && reml_priors_var.empty()) {
        const int n    = _n;
        const int k    = _woodbury_rank;
        const double Vy = _y_Ssq;          // total variance (var(y) scale used by init)
        const double trK  = _dk.sum() + static_cast<double>(n - k) * _lambda_tail;
        // tr(K²) = Σ d_j² (top-k, exact) + (n−k)·E[d_bulk²]
        // E[d_bulk²] = var(d_bulk) + mean(d_bulk)² = _tail_d_var + λ_tail²
        const double trK2 = _dk.squaredNorm()
                          + static_cast<double>(n - k) * (_tail_d_var + _lambda_tail * _lambda_tail);
        // y'Ky via woodbury_Kv: K·y = λ_tail·y + U_k·(δ⊙(U_k^T·y)), O(nk)
        const eigenVector Ky = woodbury_Kv(_y);
        const double yKy = _y.dot(Ky);
        const double yy  = _y.squaredNorm();   // y'y (≈ n·Vy when mean ≈ 0)
        const double denom = static_cast<double>(n) * trK2 - trK * trK;
        double sg_he = (denom > 1e-10)
            ? (static_cast<double>(n) * yKy - trK * yy) / denom
            : Vy * 0.5;
        // Clamp: both components must be strictly positive
        sg_he = std::max(sg_he, 0.01 * Vy);
        sg_he = std::min(sg_he, 0.99 * Vy);
        varcmp(0) = sg_he;
        varcmp(1) = std::max(Vy - sg_he, 0.01 * Vy);
    }
}

double gcta::lgL_reduce_mdl(bool no_constrain) {
    if (_r_indx.size() - 1 == 0) return 0;
    bool multi_comp = (_r_indx.size() - _r_indx_drop.size() > 1);
    LOGGER << "\nCalculating the logLikelihood for the reduced model ...\n(variance component" << (multi_comp ? "s " : " ");
    for (int i = 0; i < _r_indx.size() - 1; i++) {
        if (find(_r_indx_drop.begin(), _r_indx_drop.end(), _r_indx[i]) == _r_indx_drop.end()) LOGGER << _r_indx[i] + 1 << " ";
    }
    LOGGER << (multi_comp ? "are" : "is") << " dropped from the model)" << std::endl;
    std::vector<int> vi_buf(_r_indx);
    _r_indx = _r_indx_drop;
    eigenMatrix Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(_r_indx.size(), _r_indx.size());
    eigenVector Py(_n);
    eigenVector varcmp;
    std::vector<double> reml_priors_var, reml_priors;
    init_varcomp(reml_priors_var, reml_priors, varcmp);
    double lgL = reml_iteration(Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, false, no_constrain);
    _r_indx = vi_buf;
    return lgL;
}

double gcta::reml_iteration(eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp, bool prior_var_flag, bool no_constrain, bool reml_bivar_fix_rg)
{
    /*if(reml_bivar_fix_rg){
        if(no_constrain){
            no_constrain=false;
            LOGGER<<"Warning: --reml-no-constrain disabled. The genetic correlation is fixed so that all the variance components are constrained to be positive."<<std::endl;
        }
    }*/

    //char *mtd_str[3] = {"AI-REML", "Fisher-scoring REML", "EM-REML"};
    std::vector<std::string> mtd_str = {"AI-REML", "Fisher-scoring REML", "EM-REML"};
    int i = 0, constrain_num = 0, iter = 0, reml_mtd_tmp = _reml_mtd;
    double logdet = 0.0, logdet_Xt_Vi_X = 0.0, prev_lgL = -1e20, lgL = -1e20, dlogL = 1000.0;
    eigenVector prev_prev_varcmp(varcmp), prev_varcmp(varcmp), varcomp_init(varcmp);
    if (_reml_trace_approx && !_Vi_use_woodbury && !_bivar_reml && !_within_family) {
        LOGGER << "Using Hutch++ stochastic trace estimator with "
               << _reml_trace_approx_nprobes << " probes (memory-saving mode; "
               << "skips n*n P matrix materialisation)." << std::endl;
    }
    bool converged_flag = false;
    for (iter = 0; iter < _reml_max_iter; iter++) {
        if (reml_bivar_fix_rg) update_A(prev_varcmp);
        if (iter == 0) {
            prev_varcmp = varcomp_init;
            if (prior_var_flag){
                if(_reml_fixed_var) LOGGER << "Variance components are fixed at: " << varcmp.transpose() << std::endl;
                else LOGGER << "Prior values of variance components: " << varcmp.transpose() << std::endl;
            }
            else {
                _reml_mtd = 2;
                LOGGER << "Calculating prior values of variance components by EM-REML ..." << std::endl;
            }
        }
        if (iter == 1) {
            _reml_mtd = reml_mtd_tmp;
            LOGGER << "Running " << mtd_str[_reml_mtd] << " algorithm ..." << "\nIter.\tlogL\t";
            for (i = 0; i < _r_indx.size(); i++) LOGGER << _var_name[_r_indx[i]] << "\t";
            LOGGER << std::endl;
        }
        if (_bivar_reml) calcu_Vi_bivar(_Vi, prev_varcmp, logdet, iter); // Calculate Vi, bivariate analysis //very slow
        else if (_within_family) calcu_Vi_within_family(_Vi, prev_varcmp, logdet, iter); // within-family REML
        else {
            // Hutch++ and Woodbury both skip the n×n P materialisation (matrix-free P application).
            // Fisher scoring (mtd 1) always needs the full P for calcu_Hi every iteration.
            const bool skip_P_this_iter = (_reml_trace_approx || _Vi_use_woodbury) && _reml_mtd != 1;
            if (!calcu_Vi(_Vi, prev_varcmp, logdet, iter, skip_P_this_iter)){ // Calculate Vi
                LOGGER<<"Warning: V matrix is not positive-definite.\n";
                varcmp = prev_prev_varcmp;
                if(!calcu_Vi(_Vi, varcmp, logdet, iter)) LOGGER.e(0, "V matrix is not positive-definite.");
                if (_P.size() == 0) calcu_P(_Vi, Vi_X, Xt_Vi_X_i, _P);
                calcu_Hi(_P, Hi);
                Hi = 2 * Hi;
                break;
            }
        }
        // Hutch++ and Woodbury both skip the n×n P materialisation (matrix-free P application).
        // Fisher scoring (mtd 1) always needs the full P for calcu_Hi every iteration.
        const bool skip_P = (_reml_trace_approx || _Vi_use_woodbury) && !_bivar_reml && !_within_family && _reml_mtd != 1;
        if (skip_P) {
            logdet_Xt_Vi_X = calcu_P_nomatrix(_Vi, Vi_X, Xt_Vi_X_i);
            _P.resize(0, 0);  // release memory
        } else {
            logdet_Xt_Vi_X = calcu_P(_Vi, Vi_X, Xt_Vi_X_i, _P);
            // _Vi is not used again until the next iteration (ai_reml and calcu_tr_PA
            // operate on _P).  Releasing it here drops the working set by 1×n².
            _Vi.resize(0, 0);
        }
        if (_reml_mtd == 0) ai_reml(_P, Hi, Py, prev_varcmp, varcmp, dlogL);
        else if (_reml_mtd == 1) reml_equation(_P, Hi, Py, varcmp);
        else if (_reml_mtd == 2) em_reml(_P, Py, prev_varcmp, varcmp, dlogL);  //slow ++
        lgL = -0.5 * (logdet_Xt_Vi_X + logdet + (_y.transpose() * Py)(0, 0));

        if(_reml_force_converge && _reml_AI_not_invertible) break;
            /*{
            if(_reml_mtd != 1){
                LOGGER<<"Warning: the information matrix is not invertible. Trying to fix the problem using the Fisher-scoring approach."<<std::endl;
                _reml_mtd = 1;
                _reml_AI_not_invertible = false;
                iter--;
                continue;
            }
            else {
                LOGGER<<"Warning: the information matrix is not invertible using the Fisher-scoring approach."<<std::endl;
                break;
            }
        }*/

        // output log
        if (!no_constrain) constrain_num = constrain_varcmp(varcmp);
        if (_bivar_reml && !_bivar_no_constrain) constrain_rg(varcmp);
        if (iter > 0) {
            LOGGER << iter << "\t" << std::fixed << LOGGER.setprecision(2) << lgL << "\t";
            for (i = 0; i < _r_indx.size(); i++) LOGGER << LOGGER.setprecision(5) << varcmp[i] << "\t";

            if (constrain_num > 0) LOGGER << "(" << constrain_num << " component(s) constrained)" << std::endl;
            else LOGGER << std::endl;
        } else {
            if (!prior_var_flag) LOGGER << "Updated prior values: " << varcmp.transpose() << std::endl;
            LOGGER << "logL: " << lgL << std::endl;
            //if(_reml_max_iter==1) LOGGER<<"logL: "<<lgL<<std::endl;
        }
        if(_reml_fixed_var){
            varcmp = prev_varcmp; 
            break;
        }
        if (constrain_num * 2 > _r_indx.size()){
            if(_reml_allow_constrain_run){
                LOGGER.w(0, "more than half of the variance components are constrained.");
            }else{
                LOGGER.e(0, "analysis stopped because more than half of the variance components are constrained. The result would be unreliable.\n You may have a try of the option --reml-no-constrain.");
            }
        }
        // added by Jian Yang on 22 Oct 2014
        //if (constrain_num == _r_indx.size()) LOGGER.e(0, "analysis stopped because all variance components are constrained. You may have a try of the option --reml-no-constrain.");

        if((_reml_force_converge || _reml_no_converge) && prev_lgL > lgL){
            varcmp = prev_varcmp;
            if (_P.size() == 0) calcu_P(_Vi, Vi_X, Xt_Vi_X_i, _P);
            calcu_Hi(_P, Hi);
            Hi = 2 * Hi;
            break;
        }

        // convergence
        dlogL = lgL - prev_lgL;
        const double lgL_scale = std::max(1.0, fabs(lgL));
        const double dlogL_rel = fabs(dlogL) / lgL_scale;
        const double vc_rel = (varcmp - prev_varcmp).squaredNorm() / std::max(1.0, varcmp.squaredNorm());
        const bool like_small = (fabs(dlogL) < 1e-4)
            || (dlogL_rel < 1e-6);
        if (vc_rel < 1e-8 && like_small) {
            converged_flag = true;
            if (_reml_mtd == 2) {
                if (_P.size() == 0) calcu_P(_Vi, Vi_X, Xt_Vi_X_i, _P);
                calcu_Hi(_P, Hi);
                Hi = 2 * Hi;
            } // for calculation of SE
            break;
        }
        // Park _P's buffer into _Vi so that:
        //  (1) _P is logically empty (size 0) when calcu_Vi runs next iteration,
        //      avoiding the simultaneous _A + _Vi + _P = 3×n² ≈ 5.77 GB peak.
        //  (2) _Vi holds an n×n allocation, so calcu_Vi's resize(_n,_n) is a
        //      no-op and fills in-place — zero malloc/free overhead per iteration.
        // All code paths that use _P after ai/em_reml (calcu_Hi in break
        // branches above) exit via break before reaching here, so _P is valid
        // in those paths.  The existing _P.resize(0,0) after the loop handles
        // final cleanup.
        _Vi.swap(_P);
        prev_prev_varcmp = prev_varcmp;
        prev_varcmp = varcmp;
        prev_lgL = lgL;
    }
    
    if(_reml_fixed_var) LOGGER.w(0, "the model is evaluated at fixed variance components. The (log-)likelihood might not be maximised.");
    else {
        if(converged_flag) LOGGER << "Log-likelihood ratio converged." << std::endl;
        else {
            if(_reml_force_converge || _reml_no_converge) LOGGER.w(0, "Log-likelihood not converged. Results are not reliable.");
            else if(iter == _reml_max_iter){
                std::stringstream errmsg;
                errmsg << "Log-likelihood not converged (stop after " << _reml_max_iter << " iteractions). \nYou can specify the option --reml-maxit to allow for more iterations." << std::endl;
                if (_reml_max_iter > 1) LOGGER.e(0, errmsg.str());
            }
        }
    }
    // _P is not needed after REML convergence; free it now so it does not
    // linger through the MLMA phase (~n² doubles = ~2 GB for n≈15 500).
    _P.resize(0, 0);
    return lgL;
}

// Helper: sum all entries of Hi over an index set (replaces manual double loops)
static double sumSubmatrix(const eigenMatrix &Hi, const std::vector<int> &pos) {
    double total = 0.0;
    for (int i : pos)
        total += Hi.row(i)(Eigen::Map<const Eigen::VectorXi>(pos.data(), pos.size())).sum();
    return total;
}

// Or with Eigen 3.4+ slicing (cleaner):
// return Hi(pos, pos).sum();  // using Eigen::indexing


void gcta::calcu_Vp(double &Vp, double &Vp2, double &VarVp, double &VarVp2,
                    const eigenVector &varcmp, const eigenMatrix &Hi) {
    Vp = Vp2 = VarVp = VarVp2 = 0.0;

    auto accumulate = [&](const std::vector<int> &pos, double &vp, double &varVp) {
        for (int ii : pos) {
            vp += varcmp[ii];
            for (int jj : pos)
                varVp += Hi(ii, jj);
        }
    };

    if (_bivar_reml) {
        accumulate(_bivar_pos[0], Vp,  VarVp);
        accumulate(_bivar_pos[1], Vp2, VarVp2);
        return;
    }

    // Non-bivariate: treat _r_indx as the index set
    // varcmp[0..n-1] and full Hi block
    const Eigen::Index n = static_cast<Eigen::Index>(_r_indx.size());
    Vp    = varcmp.head(n).sum();
    VarVp = Hi.topLeftCorner(n, n).sum();
}


void gcta::calcu_hsq(int i, double Vp, double Vp2, double VarVp, double VarVp2,
                     double &hsq, double &var_hsq,
                     const eigenVector &varcmp, const eigenMatrix &Hi) {
    const double V1    = varcmp[i];
    const double VarV1 = Hi(i, i);

    // Compute hsq and var_hsq given a resolved Vp and covariance row sum
    auto compute_hsq = [&](double vp, double varVp, double cov12) {
        const double ratio = V1 / vp;
        hsq     = ratio;
        var_hsq = ratio * ratio * (VarV1 / (V1 * V1)
                                 + varVp / (vp * vp)
                                 - (2.0 * cov12) / (V1 * vp));
    };

    if (_bivar_reml) {
        // Try each bivariate partition
        auto try_partition = [&](const std::vector<int> &pos, double vp, double varVp) -> bool {
            auto iter = std::find(pos.begin(), pos.end(), i);
            if (iter == pos.end()) return false;

            // Sum Hi row i over the partition indices (Eigen row slice)
            double cov12 = 0.0;
            for (int j : pos) cov12 += Hi(*iter, j);
            // Eigen 3.4+: cov12 = Hi.row(*iter)(pos).sum();

            compute_hsq(vp, varVp, cov12);
            return true;
        };

        if (try_partition(_bivar_pos[0], Vp,  VarVp) ||
            try_partition(_bivar_pos[1], Vp2, VarVp2))
            return;

        hsq = var_hsq = -2.0;
        return;
    }

    // Non-bivariate: covariance is sum of row i across all components
    const double cov12 = Hi.row(i).head(static_cast<Eigen::Index>(_r_indx.size())).sum();
    compute_hsq(Vp, VarVp, cov12);
}


void gcta::calcu_sum_hsq(double Vp, double VarVp,
                         double &sum_hsq, double &var_sum_hsq,
                         const eigenVector &varcmp, const eigenMatrix &Hi) {
    // All components except the last (residual/error)
    const Eigen::Index n = static_cast<Eigen::Index>(_r_indx.size()) - 1;

    const double V1    = varcmp.head(n).sum();
    const double VarV1 = Hi.topLeftCorner(n, n).sum();

    // Cov12: row sums of Hi for the genetic components against ALL components
    const Eigen::Index total = static_cast<Eigen::Index>(_r_indx.size());
    const double cov12 = Hi.topLeftCorner(n, total).sum();

    const double ratio = V1 / Vp;
    sum_hsq     = ratio;
    var_sum_hsq = ratio * ratio * (VarV1 / (V1 * V1)
                                 + VarVp / (Vp * Vp)
                                 - (2.0 * cov12) / (V1 * Vp));
}

// Woodbury low-rank approximation helpers.
// K * v ≈ λ_tail·v + U_k·((D_k − λ_tail) ⊙ (U_k^T·v))
eigenVector gcta::woodbury_Kv(const eigenVector &v) const {
    eigenVector Ukv = _Uk.transpose() * v;
    Ukv.array() *= (_dk.array() - static_cast<eigenVector::Scalar>(_lambda_tail));
    return static_cast<eigenVector::Scalar>(_lambda_tail) * v + _Uk * Ukv;
}

// V^{-1} v = (v − U_k·(c_k ⊙ (U_k^T·v))) / σ²_eff
eigenVector gcta::woodbury_Viv(const eigenVector &v) const {
    eigenVector Ukv = _Uk.transpose() * v;
    Ukv.array() *= _ck.array();
    return (v - _Uk * Ukv) / static_cast<eigenVector::Scalar>(_sigma2_eff);
}

// V^{-1} Z  (matrix version)
eigenMatrix gcta::woodbury_ViZ(const eigenMatrix &Z) const {
    eigenMatrix UkZ = _Uk.transpose() * Z;
    UkZ.array().colwise() *= _ck.array();
    return (Z - _Uk * UkZ) / static_cast<eigenMatrix::Scalar>(_sigma2_eff);
}

// Compute Woodbury low-rank basis from _A[_r_indx[0]] (GRM K).
//
// Fixed-k mode (buffer_factor == 0): run SVD to rank k exactly.
// Auto-k mode  (buffer_factor  > 0): run SVD to k_max (default min(n−1, 2000)),
//   count eigenvalues above the Marchenko-Pastur bulk edge
//     λ+ = (1 + √(n/M))²
//   where M = SNP count from _grm_N or _include.size(), then set
//     k = min(k_max, max(20, ⌈k_signal × buffer_factor⌉)).
//   Using buffer_factor > 1 retains extra "transition-region" eigenpairs
//   that are treated exactly rather than collapsed to λ_tail·I, which
//   directly reduces the Jensen logdet bias in the transition band.
//
// In both modes, tr(K²) is computed to derive _tail_d_var = var(d_{k+1:n}),
// enabling a second-order logdet correction in calcu_Vi.

//TODO: can explore Nyström single-pass for woodbury basis.
void gcta::compute_woodbury_basis(int k, double buffer_factor, int k_max) {
    if (_reml_mtd == 1)
        LOGGER.e(0, "--reml-woodbury is incompatible with Fisher-scoring REML (--reml-alg 1). Use AI-REML (default) or EM-REML (--reml-alg 2).");
    if ((int)_r_indx.size() != 2)
        LOGGER.e(0, "--reml-woodbury currently supports only single-GRM models (one genetic + one residual component).");
    if (_A[_r_indx[0]].size() == 0)
        LOGGER.e(0, "--reml-woodbury: GRM component appears to be an identity matrix; cannot compute low-rank basis.");

    const bool auto_k = (buffer_factor > 0.0);
    const int  n      = _n;

    // Determine how many eigenpairs to compute
    int k_svd;
    if (auto_k) {
        k_svd = (k_max > 0) ? k_max : std::min(n - 1, 1200);
        LOGGER << "\nComputing Woodbury basis (auto-k, k_max=" << k_svd
               << ", buffer=" << buffer_factor << ") ..." << std::endl;
    } else {
        k_svd = k;
        LOGGER << "\nComputing Woodbury low-rank basis (k=" << k << ") ..." << std::endl;
    }
    if (k_svd >= n) LOGGER.e(0, "--reml-woodbury rank must be < n.");

    // Alias (double-precision mode) or upcast (single-precision mode) to double.
    // In double-precision mode eigenMatrix is already MatrixXd, so we take a
    // const reference to avoid duplicating the O(n²) GRM in memory.
#ifdef SINGLE_PRECISION
    Eigen::MatrixXd K_dbl = _A[_r_indx[0]].cast<double>();
#else
    const Eigen::MatrixXd& K_dbl = _A[_r_indx[0]];
#endif

    // ---- Randomised eigendecomposition ----
    const int oversample = 20;
    const int k_ext = std::min(k_svd + oversample, n - 1);
    // Warm-start: when mlma_loco_v2 runs successive per-chromosome REMLs, the
    // eigenvectors of G_loco_c are very close to those of G_loco_{c-1} (each chr
    // contributes ~1/n_chr of total variance).  Seeding omega with the previous
    // _Uk is equivalent to one free power-iteration step and typically cuts the
    // number of DSYMM calls from 3 → 1 for chromosomes 2..n_chr.
    // Nyström uses omega differently (C = Ω^T Y is not a subspace basis), so
    // warm-starting there is not straightforward and is skipped.
    Eigen::MatrixXd omega;
    {
        const int k_prev = static_cast<int>(_Uk.cols());
        const bool has_warm = (!_woodbury_nystrom && k_prev > 0 && _Uk.rows() == n);
        if (has_warm) {
            omega.resize(n, k_ext);
            const int k_copy = std::min(k_prev, k_ext);
            omega.leftCols(k_copy)  = _Uk.leftCols(k_copy).cast<double>();
            if (k_copy < k_ext)
                omega.rightCols(k_ext - k_copy) =
                    Eigen::MatrixXd::Random(n, k_ext - k_copy);
        } else {
            omega = Eigen::MatrixXd::Random(n, k_ext);
        }
    }
    Eigen::MatrixXd Y = K_dbl.selfadjointView<Eigen::Lower>() * omega;  // 1st (and for Nyström, only) DSYMM

    Eigen::VectorXd eval_full;
    Eigen::MatrixXd evec_full;

    if (_woodbury_nystrom) {
        // ---- Nyström single-pass (1 DSYMM total) ----
        // K ≈ Y C⁻¹ Yᵀ,  C = Ωᵀ Y = Ωᵀ K Ω  (k_ext × k_ext, symmetric)
        // Let Z = Y C^{-1/2}.  Then ZZᵀ ≈ K  and  K ≈ U S² Uᵀ  via thin SVD of Z.
        // C^{-1/2} via eigensolver: C = V Λ Vᵀ  →  C^{-1/2} = V Λ^{-1/2} Vᵀ
        // Z = W Vᵀ  where  W = Y V Λ^{-1/2}.
        // Key: SVD(Z) and SVD(W) share singular values and left singular vectors
        //      because right-multiplying by orthogonal V doesn't change them.
        // So SVD(W) directly gives eigenvalues (σ²) and eigenvectors (U).
        // Using eigensolver (not Cholesky) makes this robust to GRMs with small
        // negative eigenvalues from numerical noise in genotype standardisation.
        Eigen::MatrixXd C = omega.transpose() * Y;  // k_ext × k_ext
        omega.resize(0, 0);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_C(C);
        if (es_C.info() != Eigen::Success)
            LOGGER.e(0, "Woodbury Nyström: eigendecomposition of sketch matrix C = ΩᵀKΩ failed. "
                        "Omit --reml-woodbury-nystrom to use the default power-iteration path.");
        const double lam_max = es_C.eigenvalues().maxCoeff();
        const double eps_C   = 1e-8 * std::max(lam_max, 1.0);
        Eigen::VectorXd lam_sqrt_inv =
            es_C.eigenvalues().cwiseMax(eps_C).cwiseSqrt().cwiseInverse();
        // W = Y V diag(λ^{-1/2})  — form in-place to avoid extra n×k allocation
        Y = Y * (es_C.eigenvectors() * lam_sqrt_inv.asDiagonal());
        Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU> svd(Y);
        Y.resize(0, 0);
        eval_full = svd.singularValues().head(k_svd).array().square();
        evec_full = svd.matrixU().leftCols(k_svd);
    } else {
        // ---- Power-iteration randomised SVD (Halko et al. 2011) — 5 DSYMMs total ----
        omega.resize(0, 0);  // free: not needed after initial sketch

        constexpr int power_iter = 3;
        for (int pi = 0; pi < power_iter; ++pi) {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr_pi(Y);
            Y = qr_pi.householderQ() * Eigen::MatrixXd::Identity(n, k_ext);
            Y = K_dbl.selfadjointView<Eigen::Lower>() * Y;
        }
        {
            Eigen::HouseholderQR<Eigen::MatrixXd> qr(Y);
            Y = qr.householderQ() * Eigen::MatrixXd::Identity(n, k_ext);
        }

        Eigen::MatrixXd KY = K_dbl.selfadjointView<Eigen::Lower>() * Y;
        Eigen::MatrixXd B  = Y.transpose() * KY;
        KY.resize(0, 0);

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
        if (es.info() != Eigen::Success)
            LOGGER.e(0, "Woodbury: eigendecomposition of projected GRM failed. "
                        "Try a smaller rank or use --reml-woodbury auto.");

        eval_full = es.eigenvalues().tail(k_svd).reverse();
        evec_full = Y * es.eigenvectors().rightCols(k_svd).rowwise().reverse();
    }

    // ---- Auto-k: Marchenko-Pastur upper edge ----
    if (auto_k) {
        double M = 0.0;
        if (_grm_N.rows() == n && _grm_N.cols() == n)
            M = _grm_N.diagonal().cast<double>().mean();
        else if (!_include.empty())
            M = static_cast<double>(_include.size());
        if (M <= 0.0)
            LOGGER.e(0, "--reml-woodbury auto: cannot determine SNP count for MP threshold. "
                        "Use --reml-woodbury <k> for a fixed rank.");

        const double gamma       = static_cast<double>(n) / M;
        const double lambda_plus = std::pow(1.0 + std::sqrt(gamma), 2.0);
        const int    k_signal   = static_cast<int>((eval_full.array() > lambda_plus).count());
        const int    k_buffered = std::max(20, static_cast<int>(std::ceil(k_signal * buffer_factor)));
        k = std::min(k_svd, k_buffered);
        LOGGER << "MP bulk edge λ+ = " << lambda_plus
               << " (n=" << n << ", M=" << static_cast<long long>(M) << ")"
               << ", eigenvalues above λ+ = " << k_signal
               << ", using k = " << k << " (buffer=" << buffer_factor << ")" << std::endl;
        // Warn if k_max clamps the buffer-implied rank: the score equation is sensitive
        // to λᵢ², so quadratic-form accuracy requires k ~ k_signal×√k_signal, not just
        // k_signal. Truncation means the user should increase k_max.
        if (k_buffered > k_svd) {
            LOGGER.w(0, "Woodbury auto-k: buffer-implied k=" + std::to_string(k_buffered)
                        + " exceeds k_max=" + std::to_string(k_svd) + "; clamped to k_max. "
                        "REML score accuracy degrades when k < k_signal\u00d7\u221ak_signal. "
                        "Increase k_max via: --reml-woodbury auto:"
                        + std::to_string(k_buffered) + ":" + std::to_string(buffer_factor));
        }
    }

    // Truncate to chosen k
    Eigen::VectorXd eval = eval_full.head(k);
    Eigen::MatrixXd evec = evec_full.leftCols(k);

    // ---- λ_tail from trace ----
    const double trace_K = K_dbl.diagonal().sum();
    _lambda_tail = (trace_K - eval.sum()) / static_cast<double>(n - k);
    if (_lambda_tail < 0.0) {
        LOGGER.w(0, "Woodbury: λ_tail < 0 (" + std::to_string(_lambda_tail) + "); clamped to 0. "
                    "Consider increasing rank.");
        _lambda_tail = 0.0;
    }

    // ---- Second-order logdet correction: var of tail eigenvalues ----
    // tr(K²) = Σ K_ii² + 2·Σ_{i>j} K_ij²  (reading only the lower triangle).
    // _tail_d_var = [tr(K²) − Σ_{j≤k} d_j²] / (n−k)  −  λ_tail²
    {
        double diag_sq = K_dbl.diagonal().squaredNorm();
        double off_sq  = 0.0;
        #pragma omp parallel for reduction(+:off_sq) schedule(static)
        for (int j = 0; j < n; ++j)
            off_sq += K_dbl.col(j).tail(n - j - 1).squaredNorm();
        const double trace_K2    = diag_sq + 2.0 * off_sq;
        const double tail_sum_sq = trace_K2 - eval.squaredNorm();
        _tail_d_var = std::max(0.0, tail_sum_sq / static_cast<double>(n - k)
                                     - _lambda_tail * _lambda_tail);
    }

    // Clamp top-k eigenvalues to ≥ 0 (guard against small numerical negatives)
    eval = eval.cwiseMax(0.0);

    // Store
    _dk = eval.cast<eigenVector::Scalar>();
    _Uk = evec.cast<eigenMatrix::Scalar>();
    _woodbury_rank = k;

    // Free K to recover O(n²) memory; mark woodbury as active
    _A[_r_indx[0]].resize(0, 0);
    _Vi_use_woodbury = true;

    LOGGER << "Woodbury basis: k=" << k
           << ", λ_tail=" << _lambda_tail
           << ", tail_d_var=" << _tail_d_var << std::endl;
}

bool gcta::calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter, bool factorize_only)
{
    _Vi_use_llt = false;  // reset on every entry

    // Woodbury low-rank path: no V assembly needed.
    // Update σ²_eff, c_k, logdet analytically from the precomputed eigendecomposition.
    if (_Vi_use_woodbury) {
        if (!factorize_only)
            LOGGER.e(0, "Woodbury REML (--reml-woodbury) is incompatible with --reml-pred-rand "
                        "(explicit V^{-1} not supported in this mode).");
        Vi.resize(0, 0);
        const double sg2 = static_cast<double>(prev_varcmp[0]);
        const double se2 = static_cast<double>(prev_varcmp[_r_indx.size() - 1]);
        _sigma2_eff = sg2 * _lambda_tail + se2;
        if (_sigma2_eff <= 0.0) return false;
        _sg2 = sg2;  // cache for calcu_tr_PA_woodbury
        const int k = _woodbury_rank;
        _ck.resize(k);
        logdet = static_cast<double>(_n - k) * std::log(_sigma2_eff);
        for (int j = 0; j < k; ++j) {
            // Clamp delta to 0: if d_j < lambda_tail (can happen numerically when d_j
            // is close to the tail average), treat the contribution as zero rather than
            // introducing a negative correction that would make _ck negative.
            const double delta    = std::max(0.0, static_cast<double>(_dk[j]) - _lambda_tail);
            const double sig_delta = sg2 * delta;
            _ck[j] = static_cast<eigenVector::Scalar>(sig_delta / (_sigma2_eff + sig_delta));
            logdet += std::log(_sigma2_eff + sig_delta);
        }
        // Second-order logdet correction for tail eigenvalue variance:
        //   Exact:  Σ_{j>k} log(σ²_g·d_j + σ²_e)
        //   WB:     (n−k)·log(σ²_eff)  [Jensen over-estimates by ~(σ²_g/σ²_eff)²·(n−k)·var_tail/2]
        // This term is stored once in compute_woodbury_basis and reused each REML iteration.
        if (_tail_d_var > 0.0) {
            const double r = sg2 / _sigma2_eff;
            logdet -= 0.5 * r * r * static_cast<double>(_n - k) * _tail_d_var;
        }
        return true;
    }

    // For the Hutch++ path (factorize_only=true) reuse _Vi_L's existing n×n
    // allocation by swapping it into Vi before the resize.  On the first iteration
    // _Vi_L is empty so the swap is a no-op and resize allocates normally; on
    // subsequent iterations _Vi_L holds the previous L, the swap is O(1), and
    // Vi.resize(_n,_n) is a same-dimensions no-op — eliminating the transient
    // window where both Vi and _Vi_L are simultaneously allocated (was 2×n²).
    if (factorize_only && static_cast<int>(_r_indx.size()) > 1)
        Vi.swap(_Vi_L);
    Vi.resize(_n, _n);
    if (_r_indx.size() == 1) {
        Vi.triangularView<Eigen::Lower>().setZero();
        Vi.diagonal() = eigenVector::Constant(_n, 1.0 / prev_varcmp[0]);
        logdet = _n * log(prev_varcmp[0]);
    } 
    else {

        // Accumulate only the lower triangle of V = sum_i A_i * sigma_i^2.
        // dpotrf reads only the lower triangle, and _LLT symmetrises from lower
        // after dpotri, so the upper never needs to be filled during assembly.
        const int num_comp = _r_indx.size();

        // V = sum_i sigma_i^2 * A_i  (lower triangle only; dpotrf reads only lower).
        Vi.triangularView<Eigen::Lower>().setZero();
        // Fuse the component loop inside the column loop so each column of Vi is
        // visited exactly once, accumulating all components in a single pass.
        // The original code made num_comp separate parallel passes over n columns,
        // thrashing the L2/L3 cache for large n; this layout keeps Vi.col(j) hot
        // while walking through the A_i columns in sequence.
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < _n; j++) {
            for (int ci = 0; ci < num_comp; ci++)
                if (_A[_r_indx[ci]].size() > 0)
                    Vi.col(j).tail(_n - j) += prev_varcmp[ci] * _A[_r_indx[ci]].col(j).tail(_n - j);
        }
        // Identity components (size==0) contribute only σ² to the diagonal.
        for (int ci = 0; ci < num_comp; ci++)
            if (_A[_r_indx[ci]].size() == 0)
                Vi.diagonal().array() += prev_varcmp[ci];

        // LLT-only path: factorize V but skip dpotri — used when the explicit V^{-1}
        // matrix is not needed (Hutch++ / skip_P mode).  V is assembled in the lower
        // triangle; gcta_dpotrf factors it in-place.
        //
        // Memory note: Eigen::LLT(matrix) COPIES matrix into its internal storage,
        // creating a transient 2×n² peak every iteration.  Instead we call
        // gcta_dpotrf directly on Vi.data() — factors in-place, no copy — then
        // swap Vi's heap allocation into _Vi_L (O(1) pointer exchange).  Peak is
        // now 1×n² throughout.
        if (factorize_only && !_reml_diagV_adj) {
            gcta_blas_int blas_n = static_cast<gcta_blas_int>(_n);
            if (gcta_dpotrf(blas_n, Vi.data(), blas_n) == 0) {
                // Vi's lower triangle now holds L; diagonal = log entries for logdet.
                logdet = 2.0 * Vi.diagonal().array().log().sum();
                // O(1) swap: _Vi_L takes Vi's storage, Vi takes _Vi_L's old
                // (empty) storage.  The subsequent resize(0,0) releases that slot.
                _Vi_L.swap(Vi);
                Vi.resize(0, 0);
                _Vi_use_llt = true;
                return true;
            }
            // V not positive-definite: dpotrf has partially overwritten Vi's lower
            // triangle.  Reassemble V before falling through to full inversion.
            Vi.triangularView<Eigen::Lower>().setZero();
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < _n; j++)
                for (int ci = 0; ci < num_comp; ci++)
                    if (_A[_r_indx[ci]].size() > 0)
                        Vi.col(j).tail(_n - j) +=
                            prev_varcmp[ci] * _A[_r_indx[ci]].col(j).tail(_n - j);
            for (int ci = 0; ci < num_comp; ci++)
                if (_A[_r_indx[ci]].size() == 0)
                    Vi.diagonal().array() += prev_varcmp[ci];
        }

        constexpr bool use_lu = false;
        INVmethod method_try = _reml_inv_mtd ? static_cast<INVmethod>(_reml_inv_mtd) : (use_lu ? INV_LU : INV_LLT);
            
        /*
        for(int i = 0; i < _r_indx.size(); i++){ 
            double cur_varcmp = prev_varcmp[i];
            #pragma omp parallel for
            for(int k = 0; k < _A[_r_indx[i]].; ++k){
                for(eigenSparseMat::InnerIterator it(_A[_r_indx[i]], k); it; ++it){
                    Vi(it.row(),it.col()) += it.value() * cur_varcmp;
                }
            }
        }   
        */
        int rank = 0;
        bool ret = true;

        // Hot-path optimisation for the default LLT inversion: call dpotrf/dpotri
        // directly in-place instead of via SquareMatrixInverse.
        // SquareMatrixInverse always makes an n×n A_orig backup copy so it can restore
        // the matrix if dpotri fails and fall through to LU.  For n=15k that copy is
        // ~1.8 GB allocated on every REML iteration.  Instead, on failure we reassemble
        // V from the stored variance components (O(n²) work, O(1) extra memory) — the
        // same strategy used by the factorize_only failure path above.
        //
        // Guard with !factorize_only: if the factorize_only dpotrf above already failed
        // and fell through, Vi holds a lower-tri-only reconstruction; we let
        // SquareMatrixInverse handle it (it will fail dpotrf again and fall to LU).
        // Guard with _reml_inv_mtd==0||1: user-specified LU/QR goes through the normal
        // SquareMatrixInverse path unchanged.
        if (method_try == INV_LLT && !factorize_only) {
            gcta_blas_int blas_n_f = static_cast<gcta_blas_int>(_n);
            bool llt_ok = (gcta_dpotrf(blas_n_f, Vi.data(), blas_n_f) == 0);
            if (llt_ok) {
                logdet = Vi.diagonal().array().square().log().sum();
                llt_ok = (gcta_dpotri(blas_n_f, Vi.data(), blas_n_f) == 0);
            }
            if (llt_ok) {
                Vi.triangularView<Eigen::Upper>() = Vi.transpose();
                return true;
            }
            // V is not positive-definite (common with GRMs that have negative
            // eigenvalues).  dpotrf has partially overwritten Vi's lower triangle;
            // reassemble the full symmetric V from stored components so the LU
            // cascade below receives a valid matrix.  No extra n×n allocation needed.
            Vi.triangularView<Eigen::Lower>().setZero();
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < _n; j++)
                for (int ci = 0; ci < num_comp; ci++)
                    if (_A[_r_indx[ci]].size() > 0)
                        Vi.col(j).tail(_n - j) +=
                            prev_varcmp[ci] * _A[_r_indx[ci]].col(j).tail(_n - j);
            for (int ci = 0; ci < num_comp; ci++)
                if (_A[_r_indx[ci]].size() == 0)
                    Vi.diagonal().array() += prev_varcmp[ci];
            // Symmetrise: LU (dgetrf) requires the full matrix.
            Vi.triangularView<Eigen::Upper>() =
                Vi.triangularView<Eigen::Lower>().transpose();
            method_try = INV_LU;
            // fall through to SquareMatrixInverse below for LU/QR/SVD cascade
        }

        if(!SquareMatrixInverse(Vi, logdet, rank, method_try)){
            LOGGER<<"Warning: the variance-covaraince matrix V is invertible." << std::endl;
            if(_reml_diagV_adj == 1){
                LOGGER<<"A small positive value is added to the diagonals. The results might not be reliable!"<<std::endl;
                double d_buf = Vi.diagonal().mean() * _reml_diag_mul;
                for(int j = 0; j < _n ; j++) Vi(j,j) += d_buf;
                if(!SquareMatrixInverse(Vi, logdet, rank, method_try)){
                    LOGGER << "Still can't be inverted. Try --reml-alg-inv 2 " << std::endl;
                    ret = false;  
                }
            }else if(_reml_diagV_adj == 2){
                LOGGER << "Switching to the \"bending\" approach to invert V. This method hasn't been tested extensively. The results might not be reliable!" << std::endl;
                bend_V(Vi);
            }else{
                LOGGER.e(0, "the variance-covariance matrix V is not invertible.\nWe may try --reml-alg-inv 1 to add a small constant value to the diagonals or --reml-alg-inv 2 to bend the matrix V.");
                ret = false;
            }
                
        }
            //if(method != method_try){
        return ret;


/*
        if (!comput_inverse_logdet_LDLT(Vi, logdet)) {
            LOGGER<<"Warning: the variance-covaraince matrix V is non-positive definite." << std::endl;
            if(_reml_diagV_adj == 1){
                LOGGER<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<std::endl;
                double d_buf = Vi.diagonal().mean() * _reml_diag_mul;
                for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
                if(!comput_inverse_logdet_LDLT(Vi, logdet)){
                    LOGGER << "Still can't be inverted. Try --reml-alg-inv 2 " << std::endl;
                    return false;  
                }
            }else if(_reml_diagV_adj == 2){
                LOGGER<<"Warning: the variance-covaraince matrix V is non-positive definite." << std::endl;
                LOGGER << "Switching to the \"bending\" approach to invert V. This method hasn't been tested extensively. The results might not be reliable!" << std::endl;
                bend_V(Vi);
            }else{
                LOGGER<<"Warning: the variance-covaraince matrix V is non-positive definite." << std::endl;
                LOGGER << "Try --reml-alg-inv 1 for diagonal addition value or --reml-alg-inv 2 for bending Vi" << std::endl;
                return false;
            }
                
        }
        */

/*
        if (_V_inv_mtd == 0) {
            if (!comput_inverse_logdet_LDLT(Vi, logdet)) {
                if(_reml_force_inv) {
                    LOGGER<<"Warning: the variance-covaraince matrix V is non-positive definite." << std::endl;
                    _V_inv_mtd = 1;
                    LOGGER << "\nSwitching to the \"bending\" approach to invert V. This method hasn't been tested extensively. The results might not be reliable!" << std::endl;
                }else {
                    if(_reml_no_converge){
                        LOGGER<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<std::endl;
                        double d_buf = Vi.diagonal().mean() * 0.01;
                        for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
                        if(!comput_inverse_logdet_LDLT(Vi, logdet)) return false;  
                    } 
                    else LOGGER.e(0, "the variance-covariance matrix V is not positive definite.");
                }
            }
        }
        if (_V_inv_mtd == 1) bend_V(Vi);
        */
       /*if (_V_inv_mtd == 2) {
            if(!_reml_force_converge){
                LOGGER << "Switching from Cholesky to LU decomposition approach. The results might not be reliable!" << std::endl;
                if (!comput_inverse_logdet_LU(Vi, logdet)){
                    if(_reml_no_converge){
                        LOGGER<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<std::endl;
                        double d_buf = Vi.diagonal().mean() * 0.01;
                        for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
                        if(!comput_inverse_logdet_LDLT(Vi, logdet)) return false;  
                    } 
                    else LOGGER.e(0, "the variance-covariance matrix V is not invertible using LU decomposition.");
                }
            }
            else{
                LOGGER<<"Warning: the variance-covaraince matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<std::endl;
                double d_buf = Vi.diagonal().mean() * 0.01;
                for(j = 0; j < _n ; j++) Vi(j,j) += d_buf;
                comput_inverse_logdet_LU(Vi, logdet);
            }
        }*/
    }

    return true;
}


void gcta::bend_V(eigenMatrix &Vi)
{
    Eigen::SelfAdjointEigenSolver<eigenMatrix> eigensolver(Vi);
    eigenVector eval = eigensolver.eigenvalues();
    bending_eigenval(eval);
    eval.array() = 1.0 / eval.array();
    Vi = eigensolver.eigenvectors() * eigenDiagMat(eval) * eigensolver.eigenvectors().transpose();
}


void gcta::bend_A() {
    _Vi.resize(0, 0);
    _P.resize(0, 0);
    LOGGER << "Bending the GRM(s) to be positive-definite (may take a while if there are multiple GRMs)..." << std::endl;
    int i = 0;
    for (i = 0; i < _r_indx.size() - 1; i++) {
        #ifdef SINGLE_PRECISION
        Eigen::SelfAdjointEigenSolver<eigenMatrix> eigensolver(_A[_r_indx[i]]);
        #else
        Eigen::SelfAdjointEigenSolver<eigenMatrix> eigensolver((_A[_r_indx[i]]).cast<double>());
        #endif
        eigenVector eval = eigensolver.eigenvalues();
        if (bending_eigenval(eval)) {
            (_A[_r_indx[i]]) = eigensolver.eigenvectors() * eigenDiagMat(eval) * eigensolver.eigenvectors().transpose();
            LOGGER << "Bending the " << i + 1 << "th GRM completed." << std::endl;
        }
    }
}

bool gcta::bending_eigenval(eigenVector &eval) {
    int j = 0;
    double eval_m = eval.mean();
    if (eval.minCoeff() > 1e-6) return false;
    double S = 0.0, P = 0.0;
    for (j = 0; j < eval.size(); j++) {
        if (eval[j] >= 0) continue;
        S += eval[j];
        P = -eval[j];
    }
    //double W = S * S * 100.0 + 1;
    double W = S * S / _reml_diag_mul + 1;
    for (j = 0; j < eval.size(); j++) {
        if (eval[j] >= 0) continue;
        eval[j] = P * (S - eval[j])*(S - eval[j]) / W;
    }
    eval *= eval_m / eval.mean();
    return true;
}

bool gcta::inverse_H(eigenMatrix &H)
{    
    double d_buf = 0.0;
    INVmethod method = (_reml_inv_mtd == 0) ? INV_LLT : static_cast<INVmethod>(_reml_inv_mtd);
    int rank = 0;
    if (!SquareMatrixInverse(H, d_buf, rank, method)) return false;
    /*{
        if(_reml_force_inv) {
            LOGGER<<"Warning: the information matrix is non-positive definite. Switching from Cholesky to LU decomposition approach. The results might not be reliable!"<<std::endl;
            if (!comput_inverse_logdet_LU(H, d_buf)){
                LOGGER<<"Warning: the information matrix is invertible. A small positive value is added to the diagonals. The results might not be reliable!"<<std::endl;
                int i = 0;
                d_buf = H.diagonal().mean() * 0.001;
                for(i = 0; i < H.rows(); i++) H(i,i) += d_buf;
                if (!comput_inverse_logdet_LU(H, d_buf)) return false;
            }
        }
        else return false;
    }*/
    else return true;
}

double gcta::calcu_P(eigenMatrix &Vi, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &P)
{
    return calcu_P_impl(Vi, Vi_X, Xt_Vi_X_i, &P);
}

double gcta::calcu_P_nomatrix(eigenMatrix &Vi, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i)
{
    return calcu_P_impl(Vi, Vi_X, Xt_Vi_X_i, nullptr);
}

double gcta::calcu_P_impl(eigenMatrix &Vi, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix *P)
{
    // Vi_X = V^{-1} X.
    // When _Vi_use_llt=true, _Vi_L holds L (lower triangle from dpotrf).
    // Solve L L^T Vi_X = X via two triangular solves (DTRSM, no n×n matrix needed).
    if (_Vi_use_woodbury) {
        // Cache U_k^T X once — constant across all REML iterations.
        // woodbury_ViZ(_X) would re-compute _Uk.T*_X every call (1st of 2 O(nkc) GEMMs);
        // the _Uk_Vi_X assignment below added a 3rd.  With the cache:
        //   calcu_P_impl cost: 3×O(nkc) → 1×O(nkc) + O(kc).
        if (_UkTX.rows() != _woodbury_rank || _UkTX.cols() != (int)_X.cols())
            _UkTX.noalias() = _Uk.transpose() * _X;   // O(nkc), computed ONCE
        if (_UkTy.size() != _woodbury_rank)
            _UkTy.noalias() = _Uk.transpose() * _y;   // O(nk),  computed ONCE
        // V^{-1}X = (X − U_k diag(ck) UkTX) / σ²_eff  — single O(nkc) GEMM
        const eigenMatrix ck_UkTX = _ck.asDiagonal() * _UkTX;  // k×c, O(kc)
        Vi_X = (_X - _Uk * ck_UkTX) / static_cast<eigenMatrix::Scalar>(_sigma2_eff);
        // U_k^T V^{-1}X = (UkTX − ck_UkTX) / σ²_eff  [_Uk orthonormal ⇒ _Uk.T*_Uk = I]
        // O(kc) instead of the previous O(nkc) GEMM.
        _Uk_Vi_X = (_UkTX - ck_UkTX) / static_cast<eigenMatrix::Scalar>(_sigma2_eff);
    } else if (_Vi_use_llt) {
        Vi_X = _X;
        _Vi_L.triangularView<Eigen::Lower>().solveInPlace(Vi_X);
        _Vi_L.triangularView<Eigen::Lower>().adjoint().solveInPlace(Vi_X);
    } else {
        Vi_X.noalias() = Vi.selfadjointView<Eigen::Lower>() * _X;
    }
    Xt_Vi_X_i.noalias() = _X.transpose() * Vi_X;

    double logdet_Xt_Vi_X = 0.0;
    int rank = 0;
    INVmethod method = (_reml_inv_mtd == 0) ? INV_LLT : static_cast<INVmethod>(_reml_inv_mtd);
    if(!SquareMatrixInverse(Xt_Vi_X_i, logdet_Xt_Vi_X, rank, method)) LOGGER.e(0, "\n  the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).");

    // Cache for implicit P matvecs (applyP_vec / Hutch++ / analytical trace)
    _Vi_X = Vi_X;
    _Xt_Vi_X_i = Xt_Vi_X_i;
    // _Uk_Vi_X already set analytically inside the woodbury branch above (O(kc));
    // no additional O(nkc) GEMM needed here.

    // Skip n×n P materialisation when caller passes nullptr (Hutch++ path)
    if (!P) return logdet_Xt_Vi_X;

    // P materialisation requires the explicit V^{-1} in Vi.
    // Swap _Vi_L into Vi (O(1), no copy): Vi gets L's storage, _Vi_L becomes empty.
    // Then dpotri inverts L in-place → Vi = V^{-1}.  Peak stays at 1×n² (was 2×n²).
    if (_Vi_use_woodbury) {
        // Materialise V^{-1} = (1/σ²_eff)(I − U_k diag(c_k) U_k^T) into Vi.
        Vi.resize(_n, _n);
        Vi.setIdentity();
        Vi /= static_cast<eigenMatrix::Scalar>(_sigma2_eff);
        // Rank-k downdate: Vi -= (1/σ²_eff) U_k diag(c_k) U_k^T
        // Use scaled U_k so that rankUpdate(Z, -1) gives the correct formula.
        eigenMatrix Uk_scaled = _Uk;
        for (int j = 0; j < _woodbury_rank; ++j)
            Uk_scaled.col(j) *= std::sqrt(static_cast<double>(_ck[j]) / _sigma2_eff);
        Vi.selfadjointView<Eigen::Lower>().rankUpdate(
            Uk_scaled, static_cast<eigenMatrix::Scalar>(-1));
        Vi.triangularView<Eigen::Upper>() = Vi.transpose();
        // _Vi_use_woodbury stays true so calcu_Hi can use the K approximation.
    } else if (_Vi_use_llt) {
        Vi.swap(_Vi_L);
        gcta_blas_int blas_n_p = static_cast<gcta_blas_int>(_n);
        if (gcta_dpotri(blas_n_p, Vi.data(), blas_n_p) != 0)
            LOGGER.e(0, "dpotri failed when materialising V^{-1} for P.");
        Vi.triangularView<Eigen::Upper>() = Vi.transpose();
        _Vi_use_llt = false;  // Vi (=_Vi) is now the full V^{-1}; _Vi_L is empty
    }

    // The correction Vi_X * Xt_Vi_X_i * Vi_X^T is symmetric PSD.
    // Factor Xt_Vi_X_i = L L^T (c×c Cholesky, negligible cost since c is tiny).
    // Then Vi_X * Xt_Vi_X_i * Vi_X^T = Z * Z^T where Z = Vi_X * L (n×c).
    // selfadjointView::rankUpdate maps to BLAS dsyrk, computing only the lower
    // triangle of the n×n update — halves flops vs a general dgemm.
    Eigen::LLT<eigenMatrix> llt(Xt_Vi_X_i);
    if(llt.info() == Eigen::Success) {
        eigenMatrix Z;
        Z.noalias() = Vi_X * llt.matrixL();
        // Swap Vi → P in O(1): P inherits V^{-1}'s storage, Vi becomes empty.
        // This eliminates the simultaneous Vi+P allocation (~n²×8 ≈ 1.92 GB saved).
        // Z was computed from Vi_X before the swap, so the ordering is safe.
        P->swap(Vi);
        P->selfadjointView<Eigen::Lower>().rankUpdate(Z, static_cast<eigenMatrix::Scalar>(-1));
        P->triangularView<Eigen::Upper>() = P->transpose();
    } else {
        // Rare fallback: Xt_Vi_X_i lost positive-definiteness after inversion (numerical
        // edge case). Use general n×n products; noalias() avoids a defensive copy.
        eigenMatrix W;
        W.noalias() = Vi_X * Xt_Vi_X_i;
        P->swap(Vi);
        P->noalias() -= W * Vi_X.transpose();
    }

    return logdet_Xt_Vi_X;
}

void gcta::calcu_Hi(eigenMatrix &P, eigenMatrix &Hi) {
    // P is the REML projection matrix — symmetrise once to legitimise selfadjointView
    P = (P + P.transpose()) * 0.5;

    const int m = static_cast<int>(_r_indx.size());

    // For identity components (size==0): PA[i] = P * I = P.  Storing a full copy of P
    // (n×n, O(n²) bytes) just to compute Frobenius products is wasteful.  Instead,
    // leave PA[i] empty (size 0) and handle the identity case analytically:
    //   Hi(i,j) when i is identity  →  PA[j] ⊙ P (P symmetric, so no transpose)
    //   Hi(i,i) when both identity →  ‖P‖_F²
    // This saves 1×n² = ~1.8 GB at convergence vs the old PA[last]=P copy.
    std::vector<eigenMatrix> PA(m);
    for (int i = 0; i < m; i++) {
        if (_Vi_use_woodbury && !_bivar_reml && !_within_family && i == 0) {
            // K was freed; use K_approx = λ_tail·I + U_k·(D−λ)·U_k^T
            // PA[i] = P·K_approx = λ_tail·P + (P·U_k)·diag(D−λ)·U_k^T
            eigenMatrix PUk = P.selfadjointView<Eigen::Upper>() * _Uk;
            eigenVector delta = _dk.array() - static_cast<eigenVector::Scalar>(_lambda_tail);
            PA[i] = static_cast<eigenMatrix::Scalar>(_lambda_tail) * P;
            PA[i].noalias() += PUk * delta.asDiagonal() * _Uk.transpose();
        } else if (!_bivar_reml && !_within_family && _A[_r_indx[i]].size() == 0) {
            PA[i].resize(0, 0);  // identity marker — no copy of P
        } else
            PA[i].noalias() = P.selfadjointView<Eigen::Upper>()
                    * ((_bivar_reml || _within_family) ? _Asp[_r_indx[i]] : _A[_r_indx[i]]);
    }

    for (int i = 0; i < m; i++) {
        const bool i_id = (PA[i].size() == 0);
        for (int j = 0; j <= i; j++) {
            const bool j_id = (PA[j].size() == 0);
            double val;
            if (i_id && j_id)       val = P.squaredNorm();              // tr(P²)
            else if (i_id)          val = PA[j].cwiseProduct(P).sum();  // tr(P·PA[j])
            else if (j_id)          val = PA[i].cwiseProduct(P).sum();  // tr(PA[i]·P)
            else                    val = PA[i].cwiseProduct(PA[j].transpose()).sum();
            Hi(i, j) = Hi(j, i) = val;
        }
    }

    if (!inverse_H(Hi)) {
        if (_reml_force_converge) {
            LOGGER.w(0, "the information matrix is not invertible.");
            _reml_AI_not_invertible = true;
        } else {
            LOGGER.e(0, "the information matrix is not invertible.");
        }
    }
}

// use Fisher-scoring to estimate variance component
// input P, calculate PA, H, R and varcmp

void gcta::reml_equation(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp)
{
    // Calculate Hi
    calcu_Hi(P, Hi);
    if(_reml_AI_not_invertible) return;

    // Calculate R
    // P is symmetric; selfadjointView routes to BLAS dsymv (lower tri only).
    Py.noalias() = P.selfadjointView<Eigen::Lower>() * _y;
    eigenVector R(_r_indx.size());
    {
        // Use a reusable temporary to avoid re-allocating _n-length storage m times.
        eigenVector tmp(_n);
        for (int i = 0; i < (int)_r_indx.size(); i++) {
            if (_bivar_reml || _within_family) R(i) = (Py.transpose()*(_Asp[_r_indx[i]]) * Py)(0, 0);
            else if (_A[_r_indx[i]].size() == 0) R(i) = Py.squaredNorm();
            else {
                // selfadjointView dispatches to DSYMV (level-2 BLAS), exploiting
                // symmetry to halve memory reads vs a general DGEMV on the full matrix.
                tmp.noalias() = _A[_r_indx[i]].selfadjointView<Eigen::Lower>() * Py;
                R(i) = Py.dot(tmp);
            }
        }
    }

    // Calculate variance component
    varcmp = Hi*R;
    Hi = 2 * Hi; // for calculation of SE
}

// use Fisher-scoring to estimate variance component
// input P, calculate PA, H, R and varcmp

void gcta::ai_reml(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, double dlogL)
{
    const bool use_approx     = _reml_trace_approx && !_bivar_reml && !_within_family;
    const bool woodbury_active = _Vi_use_woodbury  && !_bivar_reml && !_within_family;

    // Py = P * y
    if (use_approx || woodbury_active) Py = applyP_vec(_y);
    else Py.noalias() = P.selfadjointView<Eigen::Lower>() * _y;
    //eigenVector cvec(_n);
    eigenMatrix APy(_n, _r_indx.size());
    //LOGGER << "AI reml 1 start" << std::endl;
    for (int i = 0; i < (int)_r_indx.size(); i++) {
        if (_bivar_reml || _within_family)
            APy.col(i).noalias() = _Asp[_r_indx[i]] * Py;
        else if (_Vi_use_woodbury && i == 0)
            APy.col(i) = woodbury_Kv(Py);
        else if (_A[_r_indx[i]].size() == 0)
            APy.col(i) = Py;
        else
            APy.col(i).noalias() = _A[_r_indx[i]].selfadjointView<Eigen::Lower>() * Py;
    }

    //LOGGER << "AI reml 2 start" << std::endl;
    // Calculate Hi
    eigenVector R(_r_indx.size());
    if (use_approx || woodbury_active) {
        // R = APy' Py — DGEMV over all m components at once
        R.noalias() = APy.transpose() * Py;
        // PAPy: one DTRSM(n, m) reads L once for all m RHS instead of m DTRSV
        // calls.  Then Hi = APy' PAPy via a single m×m DGEMM.
        const eigenMatrix PAPy = applyP_mat(APy);
        Hi.noalias() = APy.transpose() * PAPy;
    } else {
        // R(i) = Py · (A_i Py) — one DGEMV batching all m dot products
        R.noalias() = APy.transpose() * Py;
        // PAPy = P * APy — one DSYMM (level-3 BLAS) instead of m separate
        // DSYMV (level-2) calls.  n×m right-hand side lets BLAS use its blocked
        // cache-oblivious kernels, giving 2–5× better throughput for small m.
        eigenMatrix PAPy(_n, static_cast<int>(_r_indx.size()));
        PAPy.noalias() = P.selfadjointView<Eigen::Lower>() * APy;
        // Hi = APy^T * PAPy = APy^T P APy (symmetric, B^T P B for SPD P).
        // One m×m DGEMM instead of m(m+1)/2 individual DOT calls.
        Hi.noalias() = APy.transpose() * PAPy;
    }
    //LOGGER << "AI reml 2 end" << std::endl;
    Hi = 0.5 * Hi;

    // Calculate tr(PA) and dL
    eigenVector tr_PA;
    // Woodbury path: exact analytical trace O(kc²) — replaces stochastic Hutch++.
    // _Uk_Vi_X is set in calcu_P_impl each iteration.
    if (_Vi_use_woodbury && !_bivar_reml && !_within_family) {
        calcu_tr_PA_woodbury(tr_PA);
    } else if (use_approx) {
        const int eff_nprobes = (std::fabs(dlogL) < 1.0)
            ? _reml_trace_approx_nprobes
            : std::max(_reml_trace_approx_nprobes / 3, 9);
        calcu_tr_PA_hutchpp(tr_PA, eff_nprobes);
    } else calcu_tr_PA(P, tr_PA);
    R = -0.5 * (tr_PA - R);

    // Calculate variance component
    if (!inverse_H(Hi)){
        if(_reml_force_converge){
            LOGGER.w(0, "the information matrix is not invertible.");
            _reml_AI_not_invertible = true;
            return;
        }
        else LOGGER.e(0, "the information matrix is not invertible.");
    }

    eigenVector delta(_r_indx.size());
    delta = Hi*R;
    if (dlogL > 1.0) varcmp = prev_varcmp + 0.316 * delta;
    else varcmp = prev_varcmp + delta;
}

// input P, calculate varcmp

void gcta::em_reml(eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, double dlogL)
{
    // Calculate trace(PA)
    //LOGGER << "Before em_reml: " << getVMemKB() << " " << getMemKB() << ", "; 
    eigenVector tr_PA;
    if (_Vi_use_woodbury && !_bivar_reml && !_within_family) {
        // Woodbury path: exact analytical trace, no stochastic probes needed.
        calcu_tr_PA_woodbury(tr_PA);
        Py = applyP_vec(_y);
    } else if (_reml_trace_approx && !_bivar_reml && !_within_family) {
        // Adaptive probe count: use 1/3 of the budget when dlogL is large (early
        // iterations where varcmp moves far anyway), full budget near convergence
        // (|dlogL|<1) for a 3× more accurate trace estimate at the final point.
        const int eff_nprobes = (std::fabs(dlogL) < 1.0)
            ? _reml_trace_approx_nprobes
            : std::max(_reml_trace_approx_nprobes / 3, 9);
        calcu_tr_PA_hutchpp(tr_PA, eff_nprobes);
        Py = applyP_vec(_y);
    } else {
        calcu_tr_PA(P, tr_PA);  // exact path
        // P is symmetric; selfadjointView routes to BLAS dsymv, reading only the
        // lower triangle — halves memory bandwidth vs a general dgemv.
        Py.noalias() = P.selfadjointView<Eigen::Lower>() * _y;
    }
    // R(i) = Py' * A_i * Py  — computed as a dsymv + dot to avoid forming the
    // n×n outer product.  Sequential over components (typically 2–3) so BLAS
    // threads inside dsymv are not competing with an outer OMP parallel region.
    eigenVector R(_r_indx.size());
    eigenVector tmp(_n);
    for (auto i = 0; i < (int)_r_indx.size(); i++) {
        if (_bivar_reml || _within_family)
            R(i) = Py.dot(_Asp[_r_indx[i]] * Py);
        else if (_Vi_use_woodbury && i == 0) {
            tmp = woodbury_Kv(Py);
            R(i) = Py.dot(tmp);
        } else if (_A[_r_indx[i]].size() == 0) R(i) = Py.squaredNorm();
        else {
            tmp.noalias() = _A[_r_indx[i]].selfadjointView<Eigen::Lower>() * Py;
            R(i) = Py.dot(tmp);
        }
        varcmp(i) = prev_varcmp(i) - prev_varcmp(i) * prev_varcmp(i) * (tr_PA(i) - R(i)) / _n;
    }

    // Calculate variance component
    //for (i = 0; i < _r_indx.size(); i++) varcmp(i) = (prev_varcmp(i) * _n - prev_varcmp(i) * prev_varcmp(i) * tr_PA(i) + prev_varcmp(i) * prev_varcmp(i) * R(i)) / _n;

    // added by Jian Yang Dec 2014
    //varcmp = (varcmp.array() - prev_varcmp.array())*2 + prev_varcmp.array();        
}

// Implicit P*v = V^{-1}v - V^{-1}X (X'V^{-1}X)^{-1} X'V^{-1}v
// Uses cached _Vi, _Vi_X, _Xt_Vi_X_i from calcu_P.
eigenVector gcta::applyP_vec(const eigenVector &v) const
{
    // w = V^{-1} v: use Cholesky triangular solves when available (skip-P path),
    // otherwise use the precomputed explicit V^{-1} (dsymv on lower triangle).
    // _Vi_L holds L from in-place dpotrf; two DTRSV calls (L w = v, then L^T w = w)
    // give V^{-1} v without ever materialising V^{-1}.
    eigenVector w;
    if (_Vi_use_woodbury) {
        w = woodbury_Viv(v);
    } else if (_Vi_use_llt) {
        w = v;
        _Vi_L.triangularView<Eigen::Lower>().solveInPlace(w);
        _Vi_L.triangularView<Eigen::Lower>().adjoint().solveInPlace(w);
    } else {
        w = eigenVector(_Vi.selfadjointView<Eigen::Lower>() * v);
    }
    // a = X' w = X' V^{-1} v   (dgemv, c-dimensional)
    eigenVector a = _Vi_X.transpose() * v;
    // b = (X' V^{-1} X)^{-1} a   (tiny c×c symmetric solve)
    eigenVector b = _Xt_Vi_X_i.selfadjointView<Eigen::Lower>() * a;
    // P v = w - V^{-1} X b
    w.noalias() -= _Vi_X * b;
    return w;
}

// Batch version of applyP_vec: apply P to every column of Z simultaneously.
// Uses DTRSM (triangular solve with k RHS) or DSYMM instead of k separate
// DTRSV/DSYMV calls.  For a Cholesky factor L of size n×n (n²/2 entries), each
// DTRSV reads L once; k calls therefore read L k times from RAM.  DTRSM uses
// blocking to process all k columns while each cache-line of L is hot, reading
// L only once — giving a k× bandwidth improvement on the dominant step.
eigenMatrix gcta::applyP_mat(const eigenMatrix &Z) const
{
    // w_j = V^{-1} z_j for all j: one DTRSM (L read once) or DSYMM
    // _Vi_L holds L from in-place dpotrf; DTRSM(L) + DTRSM(L^T) give V^{-1} Z
    // without allocating a second n×n matrix.
    eigenMatrix W;
    if (_Vi_use_woodbury) {
        W = woodbury_ViZ(Z);
    } else if (_Vi_use_llt) {
        W = Z;
        _Vi_L.triangularView<Eigen::Lower>().solveInPlace(W);
        _Vi_L.triangularView<Eigen::Lower>().adjoint().solveInPlace(W);
    } else {
        W = eigenMatrix(_Vi.selfadjointView<Eigen::Lower>() * Z);
    }
    // a = X'V^{-1}Z = (V^{-1}X)'Z  (c×k tiny DGEMM)
    const eigenMatrix A = _Vi_X.transpose() * Z;
    // PZ = V^{-1}Z - (V^{-1}X)(X'V^{-1}X)^{-1}(X'V^{-1}Z)
    W.noalias() -= _Vi_X * (_Xt_Vi_X_i.selfadjointView<Eigen::Lower>() * A);
    return W;
}

// Hutch++ stochastic trace estimation: tr(P * A_i) for all variance components.
// NOTE: This is a MEMORY optimisation, not a speed optimisation. Each call does
// 3 DSYMM + 3 DTRSM(n,k) per component vs one O(n²) sweep for exact trace. The
// benefit is avoiding materialisation of the n×n P matrix, saving 8n² bytes of RAM.
// Probe vectors S and G are generated once (on first call) and cached in
// _hutchpp_S / _hutchpp_G so that the trace estimate is a smooth deterministic
// function of the variance components, allowing AI-REML to converge.
// Reference: Meyer, Musco, Musco, Woodruff (2021), "Hutch++", SIAM SOSA.

//TODO: we can add more probes when delta logL is small?
void gcta::calcu_tr_PA_hutchpp(eigenVector &tr_PA, int m_probes)
{
    const int ncomp = _r_indx.size();
    tr_PA.resize(ncomp);
    // Split budget: k probes each for range-sketch, low-rank, residual Hutchinson
    const int k = std::max(m_probes / 3, 3);

    // Lazily generate probe vectors once per REML run; reuse across iterations
    // so the trace estimate is a deterministic function of varcmp → smooth
    // optimisation landscape → Newton convergence.
    //
    // Use a dedicated RNG seeded from the problem dimensions rather than the
    // global StatFunc::rng() (which is seeded from hardware entropy and advances
    // non-deterministically as other code runs).  This makes probe vectors
    // reproducible across runs for the same dataset, eliminating run-to-run
    // variation in the stochastic trace estimate and therefore in the final
    // variance component estimates.  The seed mixes n and ncomp so different
    // problem sizes get different random directions.
    if (_hutchpp_S.rows() != _n || _hutchpp_S.cols() != k) {
        const uint32_t seed = static_cast<uint32_t>(_n) * 2654435761u
                              ^ (static_cast<uint32_t>(ncomp) * 2246822519u)
                              ^ 0x9e3779b9u;  // Knuth multiplicative hash
        std::mt19937 probe_rng(seed);
        std::uniform_int_distribution<int> coin(0, 1);
        _hutchpp_S.resize(_n, k);
        _hutchpp_G.resize(_n, k);
        for (int j = 0; j < k; j++)
            for (int r = 0; r < _n; r++) {
                _hutchpp_S(r, j) = coin(probe_rng) ? 1.0 : -1.0;
                _hutchpp_G(r, j) = coin(probe_rng) ? 1.0 : -1.0;
            }
    }

    for (int ci = 0; ci < ncomp; ci++) {
        const bool is_I = (_A[_r_indx[ci]].size() == 0);
        const auto &Ai  = _A[_r_indx[ci]];  // only read when !is_I

        // Phase A: range sketch K = (PA)^{1+q} * S, with q = _reml_trace_power_iter.
        // Each power iteration maps K → PA*K, squaring the eigenvalue ratio between
        // captured and missed directions.  Intermediate QR prevents numerical blow-up
        // (Halko, Martinsson, Tropp 2011).  For a flat bulk spectrum (typical human
        // GRM after PC correction) q=0 is usually sufficient.
        auto applyPA_mat = [&](const eigenMatrix &Z) -> eigenMatrix {
            if (_Vi_use_woodbury && ci == 0) {
                // K was freed; use low-rank approximation: K*Z = λ_tail·Z + U_k·((D−λ)⊙(U_k^T·Z))
                eigenMatrix UkZ = _Uk.transpose() * Z;
                UkZ.array().colwise() *= (_dk.array() - static_cast<eigenVector::Scalar>(_lambda_tail));
                eigenMatrix KZ = static_cast<eigenMatrix::Scalar>(_lambda_tail) * Z + _Uk * UkZ;
                return applyP_mat(KZ);
            }
            return is_I ? applyP_mat(Z)
                        : applyP_mat(eigenMatrix(Ai.selfadjointView<Eigen::Lower>() * Z));
        };
        eigenMatrix K = applyPA_mat(_hutchpp_S);
        for (int pw = 0; pw < _reml_trace_power_iter; pw++) {
            // Re-orthogonalize to maintain numerical stability before next matvec.
            Eigen::HouseholderQR<eigenMatrix> qr_pw(K);
            K = qr_pw.householderQ() * eigenMatrix::Identity(_n, k);
            K = applyPA_mat(K);
        }

        // QR factorisation of K → orthonormal basis Q (n × k)
        Eigen::HouseholderQR<eigenMatrix> qr(K);
        eigenMatrix Q = qr.householderQ() * eigenMatrix::Identity(_n, k);

        // Phase B: low-rank trace contribution.
        const eigenMatrix MQ = applyPA_mat(Q);
        const double t_lr = Q.cwiseProduct(MQ).sum();

        // Phase C: Hutchinson residual.
        const eigenMatrix MG = applyPA_mat(_hutchpp_G);

        // Project G orthogonal to Q: R = (I - QQ')G, MR = (I - QQ')MG
        const eigenMatrix QtG = Q.transpose() * _hutchpp_G;   // k×k
        const eigenMatrix R_  = _hutchpp_G - Q * QtG;         // n×k
        const eigenMatrix MR  = MG          - MQ * QtG;       // n×k

        tr_PA(ci) = t_lr + R_.cwiseProduct(MR).sum() / k;
    }
}

// input P, calculate tr(PA)
void gcta::calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA) {
    //LOGGER << "Before calcu_tr_PA: " << getVMemKB() << " " << getMemKB() << ", "; 

    // Calculate trace(PA)
    tr_PA.resize(_r_indx.size());
    for (int i = 0; i < _r_indx.size(); i++) {
        //LOGGER << "calcu_tr_PA " << i << std::endl;
        if (_bivar_reml || _within_family){
            //eigenMatrix temp = P * (_Asp[_r_indx[i]]);
            //LOGGER << "   matrix product finished" << std::endl;
            //tr_PA(i) = (P * (_Asp[_r_indx[i]])).diagonal().sum();   ///extremely slow
            int cur_r_indx_size = _Asp[_r_indx[i]].outerSize();
            Eigen::VectorXd v(cur_r_indx_size);
            v.setZero(cur_r_indx_size);
            #pragma omp parallel for
            for(int k = 0; k < cur_r_indx_size; ++k){
                for(eigenSparseMat::InnerIterator it(_Asp[_r_indx[i]], k); it; ++it){
                    v(k) += P(it.col(),it.row()) * it.value();
                }
            }
            tr_PA(i) = v.sum();
            v.resize(0);
            //tr_PA(i) = temp.diagonal().sum();   ///extremely slow
            //LOGGER << "diag sum finished" << std::endl;
            //temp.resize(0,0);
        }
        else if (_A[_r_indx[i]].size() == 0) {
            // Identity component: tr(P * I) = tr(P) — only the diagonal is needed.
            double s = 0.0;
            #pragma omp parallel for reduction(+:s) schedule(static)
            for (int col = 0; col < _n; col++) s += P(col, col);
            tr_PA(i) = s;
        } else {
            // tr(P * A) = P ⊙ A elementwise summed. Both P and A are symmetric,
            // so only the lower triangle is needed: diagonal contributes weight 1,
            // off-diagonal weight 2. Column-major loop order (col outer, row inner)
            // keeps both P and A accesses sequential in memory, halving bandwidth
            // vs reading all n² elements with a full ddot.
            const auto &Ai = _A[_r_indx[i]];
            double s = 0.0;
            // Replace the scalar inner loop with a BLAS DDOT via Eigen's .dot().
            // Eigen dispatches .dot() to a vectorised (SIMD) reduction, giving the
            // same arithmetic as the hand-written loop at potentially 4–8× the
            // throughput.  The outer OMP reduction is preserved for multi-threaded
            // parallelism across columns.
            // Column col has inner-loop length (n-col-1): col 0 is n-1 long, last is 0.
            // schedule(static) gives thread 0 roughly twice the work of the last thread.
            // schedule(guided) rebalances by assigning chunks of decreasing size.
            #pragma omp parallel for reduction(+:s) schedule(guided)
            for (int col = 0; col < _n; col++) {
                const int tail = _n - col - 1;
                s += P(col, col) * Ai(col, col);
                if (tail > 0)
                    s += 2.0 * P.col(col).tail(tail).dot(Ai.col(col).tail(tail));
            }
            tr_PA(i) = s;
        }
    }
}

// Exact analytical tr(PA_i) for the Woodbury low-rank path.
//
// Since V = σ²_eff·I + U_k·diag(σ²_g·δ_k)·U_k^T and
//       K_approx = λ_tail·I + U_k·diag(δ_k)·U_k^T,
// the traces can be computed algebraically without stochastic probes:
//
//   tr(V^{-1}·K_approx) = n·λ/σ²_eff + (σ²_e/(σ²_g·σ²_eff))·Σ c_k
//   tr(V^{-1})          = (n − Σ c_k) / σ²_eff
//
// The projection correction uses the cached _Vi_X (= V^{-1}X, n×c) and
// _Uk_Vi_X (= U_k^T·V^{-1}X, k×c) so no additional O(nk) GEMMs are needed:
//
//   tr(PA_i) = tr(V^{-1}·A_i) − tr(_Xt_Vi_X_i · V_X^T·A_i·V_X)
//
// Complexity: O(k·c²) after _Uk_Vi_X is set in calcu_P_impl.  With k ≈ 300
// and c ≈ 10 covariates this is negligible (<1 ms per REML iteration).
void gcta::calcu_tr_PA_woodbury(eigenVector &tr_PA) {
    const int ncomp = static_cast<int>(_r_indx.size());
    tr_PA.resize(ncomp);

    const double sigma2_eff = _sigma2_eff;
    const double sg2        = _sg2;
    const double se2        = sigma2_eff - sg2 * _lambda_tail;
    const double lambda_t   = _lambda_tail;
    const int    n          = _n;
    const int    k          = _woodbury_rank;

    // ---- Scalar quantities (O(k)) ----
    const double sum_ck = _ck.sum();

    // tr(V^{-1}) = (n − Σcₖ) / σ²_eff
    const double tr_Vinv = (n - sum_ck) / sigma2_eff;

    // tr(V^{-1}·K_approx) = n·λ/σ²_eff + (σ²_e/(σ²_g·σ²_eff))·Σcₖ
    // Guard against sg2≈0 (degenerate run); should never happen under woodbury guard.
    const double tr_Vinv_K = (sg2 > 1e-15)
        ? (n * lambda_t / sigma2_eff + (se2 / (sg2 * sigma2_eff)) * sum_ck)
        : tr_Vinv * lambda_t;  // fallback: K ≈ λI → tr(V⁻¹K) ≈ λ·tr(V⁻¹)

    // ---- Projection correction (O(k·c² + c²)) using cached _Uk_Vi_X ----
    // _Uk_Vi_X = U_k^T · V^{-1}·X  (k×c), set in calcu_P_impl.
    // _Vi_X    = V^{-1}·X           (n×c)
    // _Xt_Vi_X_i = (X^T V^{-1} X)^{-1}  (c×c)

    const int c = static_cast<int>(_Vi_X.cols());

    // V_X^T · I · V_X = V^{-1}X)^T · (V^{-1}X) = Vi_X^T Vi_X  (c×c, O(nc²))
    eigenMatrix ViXTViX(c, c);
    ViXTViX.noalias() = _Vi_X.transpose() * _Vi_X;

    // tr(M⁻¹ · ViX^T·ViX) where M = _Xt_Vi_X_i = (X^T V^{-1} X)^{-1}
    const double tr_corr_I = (_Xt_Vi_X_i * ViXTViX).trace();

    // tr(PA · I) for the residual component (last)
    tr_PA(ncomp - 1) = tr_Vinv - tr_corr_I;

    // ---- K_approx component (index 0) ----
    // V_X^T · K_approx · V_X = λ·ViXTViX + UkViX^T · diag(δ) · UkViX   (c×c, O(kc²))
    // where δ_j = d_k[j] − λ_tail.
    eigenVector delta = _dk.array() - static_cast<eigenVector::Scalar>(lambda_t);
    // Scale rows of _Uk_Vi_X by δ: result is k×c
    eigenMatrix delta_UkViX = delta.asDiagonal() * _Uk_Vi_X;
    eigenMatrix ViX_K_ViX = lambda_t * ViXTViX;
    ViX_K_ViX.noalias() += _Uk_Vi_X.transpose() * delta_UkViX;  // c×c += (c×k)(k×c)

    const double tr_corr_K = (_Xt_Vi_X_i * ViX_K_ViX).trace();

    tr_PA(0) = tr_Vinv_K - tr_corr_K;
}

// blue estimate of SNP effect

void gcta::blup_snp_geno() {
    check_autosome();

    if (_mu.empty()) calcu_mu();

    int i = 0, j = 0, k = 0, col_num = _varcmp_Py.cols();
    double x = 0.0, fcount = 0.0;

    // Calcuate A matrix
    LOGGER << "Calculating the BLUP solutions to SNP effects ..." << std::endl;
    std::vector<double> var_SNP(_include.size());
    eigenMatrix b_SNP = eigenMatrix::Zero(_include.size(), col_num); // variance of each SNP, 2pq
    for (j = 0; j < _include.size(); j++) {
        var_SNP[j] = _mu[_include[j]]*(1.0 - 0.5 * _mu[_include[j]]);
        if (fabs(var_SNP[j]) < 1.0e-50) var_SNP[j] = 0.0;
        else var_SNP[j] = 1.0 / var_SNP[j];
    }

    for (k = 0; k < _include.size(); k++) {
        fcount = 0.0;
        for (i = 0; i < _keep.size(); i++) {
            if (!_snp_1[_include[k]][_keep[i]] || _snp_2[_include[k]][_keep[i]]) {
                x = _snp_1[_include[k]][_keep[i]] + _snp_2[_include[k]][_keep[i]];
                x = (x - _mu[_include[k]]);
                for (j = 0; j < col_num; j++) b_SNP(k, j) += x * _varcmp_Py(i, j);
                fcount += 1.0;
            }
        }
        for (j = 0; j < col_num; j++) b_SNP(k, j) = (b_SNP(k, j) * var_SNP[k] / fcount)*((double) _keep.size() / (double) _include.size());
    }
    output_blup_snp(b_SNP);
}

void gcta::blup_snp_dosage() {
    check_autosome();

    compact_dosage_data(); // ensure _include[j]==j and _keep[i]==i
    if (_mu.empty()) calcu_mu();

    const int n = static_cast<int>(_keep.size());
    const int m = static_cast<int>(_include.size());
    const int col_num = _varcmp_Py.cols();

    // Subtract 2p column-wise, skipping missing values
    for (int j = 0; j < m; j++) {
        const float mu_j = static_cast<float>(_mu[j]);
        for (int i = 0; i < n; i++) {
            float& d = _geno_dose(i, j);
            if (d < DOSAGE_NA) d -= mu_j;
        }
    }

    LOGGER << "Calculating the BLUP solutions to SNP effects using imputed dosage scores ... " << std::endl;

    // var_SNP[j] = col(j).squaredNorm() / (n-1), treating missing as 0
    // For simplicity use the same column loop; missing values were left at DOSAGE_NA
    // so we must mask them. Replace DOSAGE_NA entries with 0 in a temp copy for algebra.
    // (Missing values were already skipped in subtraction, so they remain at DOSAGE_NA.)
    // Build a clean float matrix with missing → 0 for linear algebra.
    Eigen::MatrixXf G(n, m);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            const float d = _geno_dose(i, j);
            G(i, j) = (d < DOSAGE_NA) ? d : 0.0f;
        }
    }

    // var_SNP[j] = colwise squared norm / (n-1)
    Eigen::VectorXd var_SNP = G.colwise().squaredNorm().cast<double>() / (n - 1.0);
    for (int j = 0; j < m; j++) {
        if (std::fabs(var_SNP[j]) < 1.0e-50) var_SNP[j] = 0.0;
        else var_SNP[j] = 1.0 / var_SNP[j];
    }

    // b_SNP = G^T * _varcmp_Py  (m × col_num), then scale each row
    eigenMatrix b_SNP = (G.cast<double>().transpose() * _varcmp_Py.cast<double>());
    for (int k = 0; k < m; k++)
        b_SNP.row(k) *= var_SNP[k] / static_cast<double>(m);

    output_blup_snp(b_SNP);
}

void gcta::output_blup_snp(eigenMatrix &b_SNP) {
    std::string o_b_snp_file = _out + ".snp.blp";
    std::ofstream o_b_snp(o_b_snp_file.c_str());
    if (!o_b_snp) LOGGER.e(0, "cannot open the file " + o_b_snp_file + " to write.");
    int i = 0, j = 0, col_num = b_SNP.cols();
    LOGGER << "Writing BLUP solutions of SNP effects for " << _include.size() << " SNPs to [" + o_b_snp_file + "]." << std::endl;
    for (i = 0; i < _include.size(); i++) {
        o_b_snp << _snp_name[_include[i]] << "\t" << _allele_ref[_include[i]] << "\t";
        for (j = 0; j < col_num; j++) o_b_snp << b_SNP(i, j) << "\t";
        o_b_snp << std::endl;
    }
    o_b_snp.close();
    LOGGER << "BLUP solutions of SNP effects for " << _include.size() << " SNPs have been saved in the file [" + o_b_snp_file + "]." << std::endl;
}

void gcta::HE_reg(std::string grm_file, bool m_grm_flag, std::string phen_file, std::string keep_indi_file, std::string remove_indi_file, int mphen) {
    // a memory-efficient HE regression that can fit multiple GRMs
    
    int i=0, j=0, k=0, l=0, r=0, c=0, ii=0, jj=0;
    std::stringstream errmsg;
    std::vector<std::string> phen_ID, grm_id, grm_files;
    std::vector< std::vector<std::string> > phen_buf; // save individuals by column
    _id_map.clear();
    
    // find out how many GRM components
    if (m_grm_flag) {
        read_grm_filenames(grm_file, grm_files, false);
    } else {
        grm_files.push_back(grm_file);
    }
    
    // number of model terms
    unsigned n_grm = grm_files.size();
    unsigned n_term = n_grm + 1; // plus intercept
    
    // Find common individuals in GRM and phenotype files
    // first read in grm.id, which determins the order of model equations
    std::vector<std::ifstream*> A_bin;
    A_bin.resize(n_grm);
    int size_grm = 0;
    for (i = 0; i < n_grm; i++) {
        if (i==0) {
            size_grm = read_grm_id(grm_files[i], grm_id, true, true);
        } else {
            int n = read_grm_id(grm_files[i], grm_id, true, true);
            if (n != size_grm) {
                LOGGER.e(0, "file [" + grm_files[i] + "] contains a different number of individuals from other GRM files.");
            }
        }
        std::string grm_binfile = grm_files[i] + ".grm.bin";
        A_bin[i] = new std::ifstream(grm_binfile.c_str(), std::ios::in | std::ios::binary);
        if ((*A_bin[i]).bad()) LOGGER.e(0, "cannot open the file [" + grm_binfile + "] to read.");
    }
    update_id_map_kp(grm_id, _id_map, _keep);

    // read phenotypes
    read_phen(phen_file, phen_ID, phen_buf, mphen);  // ignore individuals with missing phenotypes

    update_id_map_kp(phen_ID, _id_map, _keep);
    if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if (!remove_indi_file.empty()) remove_indi(remove_indi_file);
    
    // find out matched unique ID
    // model equations (yij and Aij) will be build based on the order of this unique ID std::vector uni_id, which is in the same order of grm_id
    std::vector<std::string> uni_id;
    std::map<std::string, int> uni_id_map;
    make_uni_id(uni_id, uni_id_map);
    _n = _keep.size();
    if (_n < 1) LOGGER.e(0, "no individual is in common among the input files.");
    LOGGER << _n << " individuals are in common in these files." << std::endl;
    
    // fill phenotypes to _y std::vector based on the order of uni_id
    _y.setZero(_n);
    for (i = 0; i < phen_ID.size(); i++) {
        auto iter = uni_id_map.find(phen_ID[i]);
        if (iter == uni_id_map.end()) continue;
        _y[iter->second] = atof(phen_buf[i][mphen - 1].c_str());
    }
    
    LOGGER << "\nPerforming Haseman-Elston regression ...\n" << std::endl;

    // normalise phenotype
    LOGGER << "Standardising the phenotype ..." << std::endl;
    _y.array() -= _y.mean();
    _y.array() /= sqrt(_y.squaredNorm() / (_n - 1.0));
    
    // grm_kp contains the rows to keep in order of uni_id, which is a subset of and in the same order of grm_id
    std::vector<int> grm_kp;
    StrFunc::match(uni_id, grm_id, grm_kp);
 
    
    // initialize OLS normal equations for each individual
    // to use jackknife to estimate SE of heritability estimate
    int long n_obs = 0.5*long(_n)*(long(_n)-1);
    double z_cp=0, z_sd=0, totalSS_cp=0, totalSS_sd=0;
    
    eigenMatrix Lhs;    // X'X
    eigenVector Rhs_cp; // X'(yi*yj)
    eigenVector Rhs_sd; // X'(yi-yj)^2
    Lhs.setZero(n_term, n_term);
    Rhs_cp.setZero(n_term);
    Rhs_sd.setZero(n_term);
    
    std::vector<eigenMatrix> LhsVec(_n);
    std::vector<eigenVector> RhsCpVec(_n);
    std::vector<eigenVector> RhsSdVec(_n);
    
    Lhs(0,0) = n_obs; // X'X for intercept
    for (i = 0; i < _n; ++i) {
        LhsVec[i].setZero(n_term, n_term);
        LhsVec[i](0,0) = _n-1;
        RhsCpVec[i].setZero(n_term);
        RhsSdVec[i].setZero(n_term);
        for (j = 0; j < i; ++j) {
            z_cp = _y[i]*_y[j];
            z_sd = (_y[i] - _y[j])*(_y[i] - _y[j]);
            Rhs_cp[0] += z_cp;
            Rhs_sd[0] += z_sd;
            RhsCpVec[i][0] += z_cp;
            RhsSdVec[i][0] += z_sd;
            RhsCpVec[j][0] += z_cp;
            RhsSdVec[j][0] += z_sd;
            totalSS_cp += z_cp * z_cp;
            totalSS_sd += z_sd * z_sd;
        }
    }
    
    
    // Fill GRMij into the ordinary least squares equations without reading the whole GRM(s) into memory
    LOGGER << "Constructing ordinary least squares equations ..." << std::endl;
    eigenVector aij(n_grm);
    int size = sizeof (float);
    float f_buf = 0.0;
    float grm_cp, r_cp, r_sd;
    for (i = 0, ii = 0; i < size_grm; i++) {
        if (i != grm_kp[ii]) {
            for (j = 0; j <= i; j++) {
                for (k = 0; k < n_grm; k++) {
                    (*A_bin[k]).read((char*) &f_buf, size);
                }
            }
        }
        else {
            for (j = 0, jj = 0; j <= i; j++) {
                if (j != grm_kp[jj] || j==i) {
                    for (k = 0; k < n_grm; k++) {
                        (*A_bin[k]).read((char*) &f_buf, size);
                    }
                }
                else {
                    for (k = 0; k < n_grm; k++) {
                        (*A_bin[k]).read((char*) &f_buf, size);
                        aij[k] = f_buf;
                        r = k + 1;
                        Lhs(0,r) = Lhs(r,0) += aij[k];
                        LhsVec[ii](0,r) = LhsVec[ii](r,0) += aij[k];
                        LhsVec[jj](0,r) = LhsVec[jj](r,0) += aij[k];
                        for (l = 0; l <= k; l++) {
                            c = l + 1;
                            grm_cp = aij[k] * aij[l];
                            Lhs(c,r) = Lhs(r,c) += grm_cp;
                            LhsVec[ii](c,r) = LhsVec[ii](r,c) += grm_cp;
                            LhsVec[jj](c,r) = LhsVec[jj](r,c) += grm_cp;
                        }
                        r_cp = aij[k] * _y[ii] * _y[jj];
                        r_sd = aij[k] *(_y[ii] - _y[jj])*(_y[ii] - _y[jj]);
                        Rhs_cp[r] += r_cp;
                        Rhs_sd[r] += r_sd;
                        RhsCpVec[ii][r] += r_cp;
                        RhsCpVec[jj][r] += r_cp;
                        RhsSdVec[ii][r] += r_sd;
                        RhsSdVec[jj][r] += r_sd;
                    }
                    ++jj;
                }
            }
            ++ii;
        }
    }

    for (k = 0; k < n_grm; k++) {
        (*A_bin[k]).close();
    }
    
    // compute OLS SE and p-value
    eigenMatrix invLhs = Lhs.inverse();
    eigenVector beta_cp = invLhs * Rhs_cp;
    eigenVector beta_sd = invLhs * Rhs_sd;
    
    double sse_cp  = totalSS_cp - beta_cp.dot(Rhs_cp);
    double sse_sd  = totalSS_sd - beta_sd.dot(Rhs_sd);
    long int df = n_obs - n_term;
    double vare_cp = sse_cp/df;
    double vare_sd = sse_sd/df;
    eigenVector se_cp = (invLhs.diagonal() * vare_cp).array().sqrt();
    eigenVector se_sd = (invLhs.diagonal() * vare_sd).array().sqrt();
    
    LOGGER << "\nLeft-hand side of OLS equations (X'X)\n" << Lhs  << std::endl << std::endl;
    //LOGGER << "vare_cp " << vare_cp << std::endl;
    //LOGGER << "vare_sd " << vare_sd << std::endl << std::endl;
    
    eigenVector pval_cp(n_term);
    eigenVector pval_sd(n_term);
    for (i = 0; i < n_term; ++i) {
        float t_cp=0, t_sd=0;
        if (se_cp[i] > 0.0) t_cp = fabs(beta_cp[i] / se_cp[i]);
        if (se_sd[i] > 0.0) t_sd = fabs(beta_sd[i] / se_sd[i]);
        pval_cp[i] = StatFunc::t_prob(df, t_cp, true);
        pval_sd[i] = StatFunc::t_prob(df, t_sd, true);
    }
    
    eigenVector kvec;
    kvec.setOnes(n_term);
    kvec[0] = 0;
    double beta_sum_cp = kvec.dot(beta_cp);
    double beta_sum_sd = kvec.dot(beta_sd);
    double se_sum_cp = sqrt((kvec.transpose()*invLhs*kvec * vare_cp)(0,0));
    double se_sum_sd = sqrt((kvec.transpose()*invLhs*kvec * vare_sd)(0,0));
    double pval_sum_cp = StatFunc::t_prob(df, abs(beta_sum_cp/se_sum_cp), true);
    double pval_sum_sd = StatFunc::t_prob(df, abs(beta_sum_sd/se_sum_sd), true);
    
    
    // compute jackknife SE and p-value
    eigenMatrix betaCpMat(_n, n_term);
    eigenMatrix betaSdMat(_n, n_term);
    for (i=0; i<_n; ++i) {
        eigenMatrix invLhsi = (Lhs - LhsVec[i]).inverse();
        betaCpMat.row(i) = invLhsi*(Rhs_cp - RhsCpVec[i]);
        betaSdMat.row(i) = invLhsi*(Rhs_sd - RhsSdVec[i]);
    }
    
    eigenVector ones = eigenVector::Ones(_n);
    eigenVector jk_mean_cp = betaCpMat.colwise().mean();
    eigenVector jk_mean_sd = betaSdMat.colwise().mean();
    eigenVector jk_se_cp = (betaCpMat - ones*jk_mean_cp.transpose()).colwise().squaredNorm();
    eigenVector jk_se_sd = (betaSdMat - ones*jk_mean_sd.transpose()).colwise().squaredNorm();
    
    jk_se_cp *= (_n-1.0)/double(_n);
    jk_se_sd *= (_n-1.0)/double(_n);
    
    jk_se_cp = jk_se_cp.array().sqrt();
    jk_se_sd = jk_se_sd.array().sqrt();
    
    eigenVector betaSumCp = betaCpMat.transpose().colwise().sum();
    eigenVector betaSumSd = betaSdMat.transpose().colwise().sum();
    
    betaSumCp -= betaCpMat.col(0);  // subtract intercept
    betaSumSd -= betaSdMat.col(0);  // subtract intercept
    
    double jk_sum_mean_cp = betaSumCp.mean();
    double jk_sum_mean_sd = betaSumSd.mean();
    double jk_sum_se_cp = (betaSumCp.array() - jk_sum_mean_cp).matrix().squaredNorm();
    double jk_sum_se_sd = (betaSumSd.array() - jk_sum_mean_sd).matrix().squaredNorm();
    
    jk_sum_se_cp = sqrt((_n-1.0)/double(_n)*jk_sum_se_cp);
    jk_sum_se_sd = sqrt((_n-1.0)/double(_n)*jk_sum_se_sd);
    
    eigenVector jk_pval_cp(n_term);
    eigenVector jk_pval_sd(n_term);
    for (i = 0; i < n_term; ++i) {
        float t_cp=0, t_sd=0;
        if (jk_se_cp[i] > 0.0) t_cp = fabs(jk_mean_cp[i] / jk_se_cp[i]);
        if (jk_se_sd[i] > 0.0) t_sd = fabs(jk_mean_sd[i] / jk_se_sd[i]);
        jk_pval_cp[i] = StatFunc::t_prob(df, t_cp, true);
        jk_pval_sd[i] = StatFunc::t_prob(df, t_sd, true);
    }
    
    double jk_pval_sum_cp = StatFunc::t_prob(df, abs(jk_sum_mean_cp/jk_sum_se_cp), true);
    double jk_pval_sum_sd = StatFunc::t_prob(df, abs(jk_sum_mean_sd/jk_sum_se_sd), true);
    
    
    // output
    std::stringstream ss;
    ss << "HE-CP\n";
    ss << std::left << std::setw(16) << "Coefficient" << std::setw(16) << "Estimate" << std::setw(16) << "SE_OLS" << std::setw(16) << "SE_Jackknife" << std::setw(16) << "P_OLS" << std::setw(16) << "P_Jackknife" << std::endl;
    for (i = 0; i < n_term; ++i) {
        if (i == 0) {
            ss << std::setw(16) << "Intercept";
        } else if (n_grm==1) {
            ss << std::setw(16) << "V(G)/Vp";
        } else {
            std::stringstream vi;
            vi << "V(G" << i << ")/Vp";
            ss << std::setw(16) << vi.str();
        }
        ss << std::setw(16) << beta_cp[i] << std::setw(16) << se_cp[i] << std::setw(16) << jk_se_cp[i] << std::setw(16) << pval_cp[i] << std::setw(16) << jk_pval_cp[i] << std::endl;
    }
    if (n_grm>1) ss << std::setw(16) << "Sum of V(G)/Vp" << std::setw(16) << beta_sum_cp << std::setw(16) << se_sum_cp << std::setw(16) << jk_sum_se_cp << std::setw(16) << pval_sum_cp << std::setw(16) << jk_pval_sum_cp << std::endl;
    ss << std::endl;
    ss << "HE-SD\n";
    ss << std::setw(16) << "Coefficient" << std::setw(16) << "Estimate" << std::setw(16) << "SE_OLS" << std::setw(16) << "SE_Jackknife" << std::setw(16) << "P_OLS" << std::setw(16) << "P_Jackknife" << std::endl;
    for (i = 0; i < n_term; ++i) {
        if (i == 0) {
            ss << std::setw(16) << "Intercept";
        } else if (n_grm==1) {
            ss << std::setw(16) << "V(G)/Vp";
        } else {
            std::stringstream vi;
            vi << "V(G" << i << ")/Vp";
            ss << std::setw(16) << vi.str();
        }
        ss << std::setw(16) << -0.5*beta_sd[i] << std::setw(16) << 0.5*se_sd[i] << std::setw(16) << 0.5*jk_se_sd[i] << std::setw(16) << pval_sd[i] << std::setw(16) << jk_pval_sd[i] << std::endl;
    }
    if (n_grm>1) ss << std::setw(16) << "Sum of V(G)/Vp" << std::setw(16) << -0.5*beta_sum_sd << std::setw(16) << 0.5*se_sum_sd << std::setw(16) << 0.5*jk_sum_se_sd << std::setw(16) << pval_sum_sd << std::setw(16) << jk_pval_sum_sd << std::endl;
    LOGGER << ss.str() << std::endl;
    std::string ofile = _out + ".HEreg";
    std::ofstream os(ofile.c_str());
    if (!os) LOGGER.e(0, "cannot open the file [" + ofile + "] to write.");
    os << ss.str() << std::endl;
    LOGGER << "Results from Haseman-Elston regression have been saved in [" + ofile + "]." << std::endl;
}

void gcta::HE_reg_bivar(std::string grm_file, bool m_grm_flag, std::string phen_file, std::string keep_indi_file, std::string remove_indi_file, int mphen, int mphen2) {
    // bivariate HE regression for two traits
    
    int i=0, j=0, k=0, l=0, r=0, c=0, ii=0, jj=0, t=0;
    std::stringstream errmsg;
    std::vector<std::string> phen_ID, grm_id, grm_files;
    std::vector< std::vector<std::string> > phen_buf; // save individuals by column
    _id_map.clear();
    
    // find out how many GRM components
    if (m_grm_flag) {
        read_grm_filenames(grm_file, grm_files, false);
    } else {
        grm_files.push_back(grm_file);
    }
    
    // number of model terms
    unsigned n_grm = grm_files.size();
    unsigned n_term = n_grm + 1; // plus intercept
    
    // Find common individuals in GRM and phenotype files
    // first read in grm.id, which determins the order of model equations
    std::vector<std::ifstream*> A_bin;
    A_bin.resize(n_grm);
    int size_grm = 0;
    for (i = 0; i < n_grm; i++) {
        if (i==0) {
            size_grm = read_grm_id(grm_files[i], grm_id, true, true);
        } else {
            int n = read_grm_id(grm_files[i], grm_id, true, true);
            if (n != size_grm) {
                LOGGER.e(0, "file [" + grm_files[i] + "] contains a different number of individuals from other GRM files.");
            }
        }
        std::string grm_binfile = grm_files[i] + ".grm.bin";
        A_bin[i] = new std::ifstream(grm_binfile.c_str(), std::ios::in | std::ios::binary);
        if ((*A_bin[i]).bad()) LOGGER.e(0, "cannot open the file [" + grm_binfile + "] to read.");
    }
    update_id_map_kp(grm_id, _id_map, _keep);
    
    // read phenotypes
    _bivar_reml = true;  // need this to read in phenotype of both traits using the function below
    read_phen(phen_file, phen_ID, phen_buf, mphen, mphen2);  // ignore individuals with missing phenotypes on both traits
    _bivar_reml = false;
    
    update_id_map_kp(phen_ID, _id_map, _keep);
    if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if (!remove_indi_file.empty()) remove_indi(remove_indi_file);
    
    // find out the matched unique ID for either trait
    // model equations (yij and Aij) will be build based on the order of the unique ID std::vector, which is in the same order of grm_id
    std::vector<std::string> uni_id;
    std::map<std::string, int> uni_id_map;
    make_uni_id(uni_id, uni_id_map);
    _n = _keep.size();
    if (_n < 1) LOGGER.e(0, "no individual is in common among the input files.");
    LOGGER << _n << " individuals are in common in these files." << std::endl;
    
    
    // find out matched unique ID for each trait
    // ans store phenotypes separately for each trait based on their own unique id order
    unsigned long n1 = 0, n2 = 0;
    eigenVector y1 = eigenVector::Zero(_n);
    eigenVector y2 = eigenVector::Zero(_n);
    
    std::vector<std::string> uni_id_tr1;
    std::vector<std::string> uni_id_tr2;
    
    std::vector<int> phen_kp;  // index of phenotyped individuals
    
    mphen--;
    mphen2--;

    StrFunc::match(uni_id, phen_ID, phen_kp);
    for (i=0; i < phen_kp.size(); i++) {
        int idx = phen_kp[i];
        if (phen_buf[idx][mphen] != "NA" && phen_buf[idx][mphen] != "-9") {
            uni_id_tr1.push_back(uni_id[i]);
            y1[n1++] = atof(phen_buf[idx][mphen].c_str());
            //LOGGER << i << " " << idx << " " << phen_buf[idx][mphen] << " " << atof(phen_buf[idx][mphen].c_str()) << std::endl;
        }
        if (phen_buf[idx][mphen2] != "NA" && phen_buf[idx][mphen2] != "-9") {
            uni_id_tr2.push_back(uni_id[i]);
            y2[n2++] = atof(phen_buf[idx][mphen2].c_str());
        }
        ++ii;
    }
    y1.conservativeResize(n1);
    y2.conservativeResize(n2);
    
    LOGGER << y1.size() << " non-missing phenotypes for trait #1 and " << y2.size() << " for trait #2" << std::endl;
    if (y1.size()==0) LOGGER.e(0, "no non-missing phenotypes for trait 1.");
    if (y2.size()==0) LOGGER.e(0, "no non-missing phenotypes for trait 2.");

    // grm_kp contains the rows of grm_id to keep in order of uni_id, which is a subset of and in the same order of grm_id
    std::vector<int> grm_kp_tr1;
    std::vector<int> grm_kp_tr2;
    
    StrFunc::match(uni_id_tr1, grm_id, grm_kp_tr1);
    StrFunc::match(uni_id_tr2, grm_id, grm_kp_tr2);
    
    
    LOGGER << "\nPerforming Haseman-Elston regression ...\n" << std::endl;
    
    // normalise phenotype
    LOGGER << "Standardising the phenotype ..." << std::endl;
    y1.array() -= y1.mean();
    y2.array() -= y2.mean();
    y1.array() /= sqrt(y1.squaredNorm() / (n1 - 1.0));
    y2.array() /= sqrt(y2.squaredNorm() / (n2 - 1.0));
    

    // initialize OLS normal equations for each individual
    // to use jackknife to estimate SE of heritability estimate
    std::vector<eigenMatrix> Lhs(3);  // X'X       for trait 1, 2 and their covariance
    std::vector<eigenVector> Rhs(3);  // X'(yi*yj) for trait 1, 2 and their covariance
    std::vector<std::vector<eigenMatrix> > LhsJk(3);  // Leave-one-individual-out (Jackknife) Lhs for the 3 var-cov componenets
    std::vector<std::vector<eigenVector> > RhsJk(3);  // Leave-one-individual-out (Jackknife) Rhs for the 3 var-cov componenets

    unsigned long minn = std::min(n1,n2);
    unsigned long maxn = std::max(n1,n2);

    std::vector<unsigned long> nObs = {n1*(n1-1)/2, n2*(n2-1)/2, n1*n2};
    std::vector<unsigned long> nJk  = {n1, n2, minn};
    
    double z = 0.0;
    std::vector<double> totalSS = {0, 0, 0};

    for (t = 0; t < 3; ++t) {
        Lhs[t].setZero(n_term, n_term);
        Lhs[t](0,0) = nObs[t];   // X'X for intercept
        Rhs[t].setZero(n_term);
        LhsJk[t].resize(nJk[t]);
        RhsJk[t].resize(nJk[t]);
    }
    
    for (i=0; i<n1; ++i) {  // for trait 1
        LhsJk[0][i].setZero(n_term, n_term);
        LhsJk[0][i](0,0) = n1-1;
        RhsJk[0][i].setZero(n_term);
        for (j=0; j<i; ++j) {
            z = y1[i]*y1[j];
            Rhs[0][0] += z;
            RhsJk[0][i][0] += z;
            RhsJk[0][j][0] += z;
            totalSS[0] += z*z;
        }
    }
    
    for (i=0; i<n2; ++i) {  // for trait 2
        LhsJk[1][i].setZero(n_term, n_term);
        LhsJk[1][i](0,0) = n2-1;
        RhsJk[1][i].setZero(n_term);
        for (j=0; j<i; ++j) {
            z = y2[i]*y2[j];
            Rhs[1][0] += z;
            RhsJk[1][i][0] += z;
            RhsJk[1][j][0] += z;
            totalSS[1] += z*z;
        }
    }
    
    eigenVector *ymin = minn == n1 ? &y1 : &y2;
    eigenVector *ymax = maxn == n1 ? &y1 : &y2;
    
    for (i=0; i<minn; ++i) {  // for covariance between trait 1 and trait 2
        LhsJk[2][i].setZero(n_term, n_term);
        LhsJk[2][i](0,0) = n1+n2-1;
        RhsJk[2][i].setZero(n_term);
    }
    for (i=0; i<minn; ++i) {
        for (j=0; j<maxn; ++j) {
            z = (*ymin)[i]*(*ymax)[j];
            Rhs[2][0] += z;
            if (j<minn) {
                RhsJk[2][i][0] += z;
                RhsJk[2][j][0] += z;
            }
            if (i==j) RhsJk[2][i][0] -= z;
            totalSS[2] += z*z;
        }
    }
    
    
    // Fill GRMij into the ordinary least squares equations without reading the whole GRM(s) into memory
    LOGGER << "Constructing ordinary least squares equations ..." << std::endl;
    eigenVector aij(n_grm);
    int size = sizeof (float);
    float f_buf = 0.0;
    float lhs, rhs;
    long t1i=0, t2i=0, t12i=0, t21i=0;
    long t1j=0, t2j=0, t12j=0, t21j=0;
    
    unsigned long count[] = {0,0,0};
    
    for (i = 0; i < size_grm; ++i) {
        t1j = t2j = t12j = t21j = 0;
        for (j = 0; j <= i; ++j) {
            for (k = 0; k < n_grm; k++) {
                (*A_bin[k]).read((char*) &f_buf, size);
                if (i == j) continue;
                
                if (i == grm_kp_tr1[t1i] && j == grm_kp_tr1[t1j]) {   // Trait 1
                    aij[k] = f_buf;
                    r = k + 1;   // first one is intercept
                    Lhs[0](0,r) = Lhs[0](r,0) += aij[k];   // symetric
                    LhsJk[0][t1i](0,r) = LhsJk[0][t1i](r,0) += aij[k];
                    LhsJk[0][t1j](0,r) = LhsJk[0][t1j](r,0) += aij[k];
                    for (l = 0; l <= k; l++) {
                        c = l + 1;
                        lhs = aij[k] * aij[l];
                        Lhs[0](c,r) = Lhs[0](r,c) += lhs;
                        LhsJk[0][t1i](c,r) = LhsJk[0][t1i](r,c) += lhs;
                        LhsJk[0][t1j](c,r) = LhsJk[0][t1j](r,c) += lhs;
                    }
                    rhs = aij[k] * y1[t1i] * y1[t1j];
                    Rhs[0][r] += rhs;
                    RhsJk[0][t1i][r] += rhs;
                    RhsJk[0][t1j][r] += rhs;
                    if (k == n_grm-1) {
                        ++t1j;
                        ++count[0];
                    }
                }
                
                if (i == grm_kp_tr2[t2i] && j == grm_kp_tr2[t2j]) {   // Trait 2
                    aij[k] = f_buf;
                    r = k + 1;   // first one is intercept
                    Lhs[1](0,r) = Lhs[1](r,0) += aij[k];   // symetric
                    LhsJk[1][t2i](0,r) = LhsJk[1][t2i](r,0) += aij[k];
                    LhsJk[1][t2j](0,r) = LhsJk[1][t2j](r,0) += aij[k];
                    for (l = 0; l <= k; l++) {
                        c = l + 1;
                        lhs = aij[k] * aij[l];
                        Lhs[1](c,r) = Lhs[1](r,c) += lhs;
                        LhsJk[1][t2i](c,r) = LhsJk[1][t2i](r,c) += lhs;
                        LhsJk[1][t2j](c,r) = LhsJk[1][t2j](r,c) += lhs;
                    }
                    rhs = aij[k] * y2[t2i] * y2[t2j];
                    Rhs[1][r] += rhs;
                    RhsJk[1][t2i][r] += rhs;
                    RhsJk[1][t2j][r] += rhs;
                    if (k == n_grm-1) {
                        ++t2j;
                        ++count[1];
                    }
                }
                
                if (i == grm_kp_tr1[t12i] && j == grm_kp_tr2[t21j]) {  // Trait 1 x Trait 2
                    aij[k] = f_buf;
                    r = k + 1;
                    Lhs[2](0,r) = Lhs[2](r,0) += aij[k];
                    if (t12i < minn && t21j < minn) {  // square matrix of size of min(n1, n2)
                        LhsJk[2][t12i](0,r) = LhsJk[2][t12i](r,0) += aij[k];
                        LhsJk[2][t21j](0,r) = LhsJk[2][t21j](r,0) += aij[k];  // crossout the row and the column
                    }
                    if (t12i == t21j) {
                        LhsJk[2][t12i](0,r) = LhsJk[2][t12i](r,0) -= aij[k];  // the cross-point is double counted
                    }
                    for (l = 0; l <= k; l++) {
                        c = l + 1;
                        lhs = aij[k] * aij[l];
                        Lhs[2](c,r) = Lhs[2](r,c) += lhs;
                        if (t12i < minn && t21j < minn) {
                            LhsJk[2][t12i](c,r) = LhsJk[2][t12i](r,c) += lhs;
                            LhsJk[2][t21j](c,r) = LhsJk[2][t21j](r,c) += lhs;
                        }
                        if (t12i == t21j) {
                            LhsJk[2][t12i](c,r) = LhsJk[2][t12i](r,c) -= lhs;
                        }
                    }
                    rhs = aij[k] * y1[t12i] * y2[t21j];
                    Rhs[2][r] += rhs;
                    if (t12i < minn && t21j < minn) {
                        RhsJk[2][t12i][r] += rhs;
                        RhsJk[2][t21j][r] += rhs;
                    }
                    if (t12i == t21j) {
                        RhsJk[2][t12i][r] -= rhs;
                    }
                    if (k == n_grm-1) {
                        ++t21j;
                        ++count[2];
                    }
                }
                else if (i == grm_kp_tr2[t21i] && j == grm_kp_tr1[t12j]) {  // Trait 2 x Trait 1
                    aij[k] = f_buf;
                    r = k + 1;
                    Lhs[2](0,r) = Lhs[2](r,0) += aij[k];
                    if (t21i < minn && t12j < minn) {
                        LhsJk[2][t12j](0,r) = LhsJk[2][t12j](r,0) += aij[k];
                        LhsJk[2][t21i](0,r) = LhsJk[2][t21i](r,0) += aij[k];
                    }
                    if (t21i == t12j) {
                        LhsJk[2][t12j](0,r) = LhsJk[2][t12j](r,0) -= aij[k];
                    }
                    for (l = 0; l <= k; l++) {
                        c = l + 1;
                        lhs = aij[k] * aij[l];
                        Lhs[2](c,r) = Lhs[2](r,c) += lhs;
                        if (t21i < minn && t12j < minn) {
                            LhsJk[2][t12j](c,r) = LhsJk[2][t12j](r,c) += lhs;
                            LhsJk[2][t21i](c,r) = LhsJk[2][t21i](r,c) += lhs;
                        }
                        if (t21i == t12j) {
                            LhsJk[2][t12j](c,r) = LhsJk[2][t12j](r,c) -= lhs;
                        }
                    }
                    rhs = aij[k] * y1[t12j] * y2[t21i];
                    Rhs[2][r] += rhs;
                    if (t21i < minn && t12j < minn) {
                        RhsJk[2][t12j][r] += rhs;
                        RhsJk[2][t21i][r] += rhs;
                    }
                    if (t21i == t12j) {
                        RhsJk[2][t12j][r] -= rhs;
                    }
                    if (k == n_grm-1) {
                        ++t12j;
                        ++count[2];
                    }
                }
            }
        }
        
        if (i == grm_kp_tr1[t1i] && t1i < n1-1) ++t1i;
        if (i == grm_kp_tr2[t2i] && t2i < n2-1) ++t2i;
        if (i == grm_kp_tr1[t12i]) ++t12i;
        if (i == grm_kp_tr2[t21i]) ++t21i;
    }
    
    
    LOGGER << "\n length of covariates:" << std::endl;
    LOGGER << "\t trait1:  " << count[0] << std::endl;
    LOGGER << "\t trait2:  " << count[1] << std::endl;
    LOGGER << "\t trait12: " << count[2] << std::endl;
    
    
    for (k = 0; k < n_grm; k++) {
        (*A_bin[k]).close();
    }
    
    // print X'X
    eigenMatrix LhsAll;
    LhsAll.setZero(3*n_term, 3*n_term);
    LhsAll.topLeftCorner(n_term, n_term) = Lhs[0];
    LhsAll.block(n_term, n_term, n_term, n_term) = Lhs[1];
    LhsAll.bottomRightCorner(n_term, n_term) = Lhs[2];
    
    LOGGER << "\nLeft-hand side of OLS equations (X'X)\n" << LhsAll  << std::endl << std::endl;
    
    // compute OLS variance-covariance matrix of estimates and p-value
    std::vector<eigenMatrix> invLhs(3);
    std::vector<eigenVector> beta(3);
    std::vector<eigenVector> se(3);
    std::vector<eigenVector> pval(3);
    std::vector<double> vare(3);
    std::vector<long unsigned> df(3);
    eigenVector kvec;
    kvec.setOnes(n_term);
    kvec[0] = 0;
    std::vector<double> betaSum(3);
    std::vector<double> seSum(3);
    std::vector<double> pvalSum(3);
    for (t=0; t<3; ++t) {
        invLhs[t] = Lhs[t].inverse();
        beta[t] = invLhs[t]*Rhs[t];
        double sse = totalSS[t] - beta[t].dot(Rhs[t]);
        df[t] = nObs[t] - n_term;
        vare[t] = sse/df[t];
        se[t] = (invLhs[t].diagonal() * vare[t]).array().sqrt();
        eigenVector tstat = (beta[t].array()/se[t].array()).abs();
        pval[t].setZero(n_term);
        for (i=0; i<n_term; ++i) {
            pval[t][i] = StatFunc::t_prob(df[t], tstat[i], true);
        }
        betaSum[t] = kvec.dot(beta[t]);
        seSum[t] = sqrt((kvec.transpose()*invLhs[t]*kvec * vare[t])(0,0));
        pvalSum[t] = StatFunc::t_prob(df[t], abs(betaSum[t]/seSum[t]), true);
    }
    
    
    // compute jackknife SE and p-value
    std::vector<eigenMatrix> betaJk(3);
    std::vector<eigenVector> seJk(3);
    std::vector<eigenVector> pvalJk(3);
    std::vector<eigenVector> betaSumJk(3);
    std::vector<double> seSumJk(3);
    std::vector<double> pvalSumJk(3);
    eigenMatrix invLhsi;
    for (t=0; t<3; ++t) {
        betaJk[t].setZero(nJk[t], n_term);
        for (i=0; i<nJk[t]; ++i) {
            invLhsi = (Lhs[t] - LhsJk[t][i]).inverse();
            betaJk[t].row(i) = invLhsi*(Rhs[t]-RhsJk[t][i]);
        }
        eigenVector ones = eigenVector::Ones(nJk[t]);
        eigenVector betaMeanJk = betaJk[t].colwise().mean();
        seJk[t] = (betaJk[t] - ones*betaMeanJk.transpose()).colwise().squaredNorm();
        seJk[t] *= (nJk[t]-1)/double(nJk[t]);
        seJk[t] = seJk[t].array().sqrt();
        
        betaSumJk[t] = betaJk[t].transpose().colwise().sum();
        betaSumJk[t] -= betaJk[t].col(0);  // subtract intercept
        double betaSumMeanJk = betaSumJk[t].mean();
        seSumJk[t] = (betaSumJk[t].array() - betaSumMeanJk).matrix().squaredNorm();
        seSumJk[t] = sqrt((nJk[t]-1)/double(nJk[t])*seSumJk[t]);
        
        eigenVector tstat = (betaMeanJk.array()/seJk[t].array()).abs();
        pvalJk[t].setZero(n_term);
        for (i = 0; i < n_term; ++i) {
            pvalJk[t][i] = StatFunc::t_prob(df[t], tstat[i], true);
        }
        pvalSumJk[t] = StatFunc::t_prob(df[t], abs(betaSumMeanJk/seSumJk[t]), true);
    }
    
    
    // genetic correlation between traits
    std::vector<double> rG(n_grm);
    std::vector<double> rG_se(n_grm);
    std::vector<double> rG_seJk(n_grm);
    double V1=0, V2=0, C=0, varV1=0, varV2=0, varC=0;
    for (i=0; i<n_grm; ++i) {
        V1 = beta[0][i+1];
        V2 = beta[1][i+1];
        C  = beta[2][i+1];
        varV1 = se[0][i+1]*se[0][i+1];
        varV2 = se[1][i+1]*se[1][i+1];
        varC  = se[2][i+1]*se[2][i+1];
        rG[i] = C/sqrt(V1*V2);
        rG_se[i] = sqrt(rG[i]*rG[i]*(varC/(C*C) + varV1/(4*V1*V1) + varV2/(4*V2*V2)));   // OLS estimate with Talor serial
        
        // Jackknife estimate of SE of rG
        eigenVector rgJk = betaJk[2].col(i+1).head(minn).array() / (betaJk[0].col(i+1).head(minn).array() * betaJk[1].col(i+1).head(minn).array()).sqrt();
        rG_seJk[i] = (rgJk.array() - rgJk.mean()).matrix().squaredNorm();
        rG_seJk[i] = sqrt((minn-1)/double(minn)*rG_seJk[i]);
    }
    V1 = betaSum[0];
    V2 = betaSum[1];
    C  = betaSum[2];
    varV1 = seSum[0]*seSum[0];
    varV2 = seSum[1]*seSum[1];
    varC  = seSum[2]*seSum[2];
    double rG_sum = C/sqrt(V1*V2);
    double rG_sum_se = sqrt(rG_sum*rG_sum*(varC/(C*C) + varV1/(4*V1*V1) + varV2/(4*V2*V2)));
    eigenVector rgJk = betaSumJk[2].head(minn).array() / (betaSumJk[0].head(minn).array() * betaSumJk[1].head(minn).array()).sqrt();
    double rG_sum_seJk = (rgJk.array() - rgJk.mean()).matrix().squaredNorm();
    rG_sum_seJk = sqrt((minn-1)/double(minn)*rG_sum_seJk);

    
    // sampling variance-covariance of estimates of variance components
    LOGGER << "\nJackknife sampling variance/covariance of the estimates of heritability:" << std::endl;
    eigenMatrix betaGJk(minn, 3*n_grm);
    j = 0;
    for (i=0; i<n_grm; ++i) {
        for (t=0; t<3; ++t) {
            betaGJk.col(j++) = betaJk[t].col(i+1);
        }
    }
    eigenMatrix centered = betaGJk.rowwise() - betaGJk.colwise().mean();
    eigenMatrix varcov = (centered.adjoint() * centered) / double(minn);
    varcov *= double(minn-1);
    
    LOGGER << varcov << std::endl << std::endl;
    
    
    // output
    std::vector<int> trid = {1, 2, 12};
    std::vector<std::string> trch = {"V", "V", "C"};
    std::stringstream ss;
    ss << "\nHE-CP\n";
    ss << std::left << std::setw(20) << "Coefficient" << std::setw(16) << "Estimate" << std::setw(16) << "SE_OLS" << std::setw(16) << "SE_Jackknife" << std::setw(16) << "P_OLS" << std::setw(16) << "P_Jackknife" << std::endl;
    for (t=0; t<3; ++t) {
        std::stringstream tmp;
        tmp << "Intercept_tr" << trid[t];
        ss  << std::setw(20) << tmp.str() << std::setw(16) << beta[t][0] << std::setw(16) << se[t][0] << std::setw(16) << seJk[t][0] << std::setw(16) << pval[t][0] << std::setw(16) << pvalJk[t][0] << std::endl;
    }
    for (i = 1; i < n_term; ++i) {
        for (t=0; t<3; ++t) {
            std::stringstream tmp;
            if (n_grm==1) tmp << trch[t] << "(G)/Vp_tr" << trid[t];
            else          tmp << trch[t] << "(G" << i << ")/Vp_tr" << trid[t];
            ss  << std::setw(20) << tmp.str() << std::setw(16) << beta[t][i] << std::setw(16) << se[t][i] << std::setw(16) << seJk[t][i] << std::setw(16) << pval[t][i] << std::setw(16) << pvalJk[t][i] << std::endl;
        }
    }
    if (n_grm>1) {
        for (t=0; t<3; ++t) {
            std::stringstream tmp;
            tmp << "Sum of " << trch[t] << "(G)/Vp_tr" << trid[t];
            ss << std::setw(20) << tmp.str() << std::setw(16) << betaSum[t] << std::setw(16) << seSum[t] << std::setw(16) << seSumJk[t] << std::setw(16) << pvalSum[t] << std::setw(16) << pvalSumJk[t] << std::endl;
        }
    }
    for (i = 0; i < n_grm; ++i) {
        std::stringstream tmp;
        if (n_grm==1) tmp << "rG";
        else          tmp << "rG" << i+1;
        ss << std::setw(20) << tmp.str() << std::setw(16) << rG[i] << std::setw(16) << rG_se[i] << std::setw(16) << rG_seJk[i] << std::endl;
    }
    if (n_grm>1) {
        ss << std::setw(20) << "Total rG" << std::setw(16) << rG_sum << std::setw(16) << rG_sum_se << std::setw(16) << rG_sum_seJk << std::endl;
    }
    ss << std::setw(20) << "N_tr1" << std::setw(16) << n1 << std::endl;
    ss << std::setw(20) << "N_tr2" << std::setw(16) << n2 << std::endl;
    ss << std::endl;
    LOGGER << ss.str() << std::endl;
    std::string ofile = _out + ".HEreg";
    std::ofstream os(ofile.c_str());
    if (!os) LOGGER.e(0, "cannot open the file [" + ofile + "] to write.");
    os << ss.str() << std::endl;
    LOGGER << "Results from Haseman-Elston regression have been saved in [" + ofile + "]." << std::endl;
}

/*   // old implementation of single component HE regression
void gcta::HE_reg(std::string grm_file, std::string phen_file, std::string keep_indi_file, std::string remove_indi_file, int mphen) {
    int i = 0, j = 0, k = 0, l = 0;
    std::stringstream errmsg;
    std::vector<std::string> phen_ID, grm_id;
    std::vector< std::vector<std::string> > phen_buf; // save individuals by column

    read_grm(grm_file, grm_id, true, false, true);
    update_id_map_kp(grm_id, _id_map, _keep);
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if (!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if (!remove_indi_file.empty()) remove_indi(remove_indi_file);

    std::vector<std::string> uni_id;
    std::map<std::string, int> uni_id_map;
    make_uni_id(uni_id, uni_id_map);
    _n = _keep.size();
    if (_n < 1) LOGGER.e(0, "no individual is in common among the input files.");

    _y.setZero(_n);
    for (i = 0; i < phen_ID.size(); i++) {
        auto iter = uni_id_map.find(phen_ID[i]);
        if (iter == uni_id_map.end()) continue;
        _y[iter->second] = atof(phen_buf[i][mphen - 1].c_str());
    }

    LOGGER << "\nPerforming Haseman-Elston regression ...\n" << std::endl;
    int n = _n * (_n - 1) / 2;
//       std::vector<bool> nomiss(_n*(_n-1));
//       for(i=0, k=0; i<_n; i++){
//           for(j=0; j<i; j++, k++){
//               if(CommFunc::FloatNotEqual(_grm(i,j),0)){
//                   nomiss[k]=true;
//                   n++;
//               }
//               else nomiss[k]=false;
//           }
//       }

    // normalise phenotype
    LOGGER << "Standardising the phenotype ..." << std::endl; 
    eigenVector y = _y.array() - _y.mean();
    y = y.array() / sqrt(y.squaredNorm() / (_n - 1.0));

    std::vector<double> y_cp(n), y_sd(n), x(n), rst;
    for (i = 0, k = 0; i < _n; i++) {
        for (j = 0; j < i; j++, k++) {
            y_sd[k] = (y[i] - y[j])*(y[i] - y[j]);
            y_cp[k] = y[i]*y[j];
            x[k] = _grm(i, j);   // Is this a bug? should i,j refer to ind in the keep list??
        }
    }
    eigenMatrix reg_sum_sd = reg(y_sd, x, rst, true);
    eigenMatrix reg_sum_cp = reg(y_cp, x, rst, true);

    std::stringstream ss;
    ss << "HE-SD" << std::endl;
    ss << "Coefficient\tEstimate\tSE\tP\n";
    ss << "Intercept\t" << reg_sum_sd.row(0) << std::endl;
    ss << "Slope\t" << reg_sum_sd.row(1) << std::endl;
    ss << "V(G)/Vp\t" << -1 * reg_sum_sd(1, 0) / reg_sum_sd(0.0) << "\t" << reg_sum_sd(1, 1) / fabs(reg_sum_sd(0.0)) << std::endl;
    ss << "\nHE-CP" << std::endl;
    ss << "Coefficient\tEstimate\tSE\tP\n";
    ss << "Intercept\t" << reg_sum_cp.row(0) << std::endl;
    ss << "Slope\t" << reg_sum_cp.row(1) << std::endl;
    ss << "V(G)/Vp\t" << reg_sum_cp(1, 0) << "\t" << reg_sum_cp(1, 1) << std::endl;

    LOGGER << ss.str() << std::endl;
    std::string ofile = _out + ".HEreg";
    std::ofstream os(ofile.c_str());
    if (!os) LOGGER.e(0, "cannot open the file [" + ofile + "] to write.");
    os << ss.str() << std::endl;
    LOGGER << "Results from Haseman-Elston regression have been saved in [" + ofile + "]." << std::endl;

} */
