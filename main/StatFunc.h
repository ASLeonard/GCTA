/*
 * Interface to the statistical functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#ifndef _STATFUNC_H
#define _STATFUNC_H

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>
#include "CommFunc.h"
#include "dcdflib.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

namespace StatFunc
{
    ////////// P-value Calculatiion Functions Start ////////////////
    // if two_tail=true, return two tail test probability;
    // otherwise, return one tail test probability;
    double t_prob(double df, double t_value, bool two_tail = true);

    double F_prob(double df_1, double df_2, double F_value);

    double betai(double a, double b, double x);

    double gammln(double xx);

    double betacf(double a, double b, double x);

    double chi_prob(double df, double chi_sqr_val);

    double gammp(const double a, const double x);

    void gser(double &gamser, const double a, const double x, double &gln);

    void gcf(double &gammcf, const double a, const double x, double &gln);
    ////////// P-value Calculatiion Functions End ////////////////

    ///////// Random Number Generation Functions Start ////////
    // (0, 1) uniform distribution generator
    double ran1(int &idum);

    // (a, b) uniform distribution generator
    double UniformDev(double a, double b, int &idum);

    // normal distribution generator
    double gasdev(int &idum);

    // generate a sequence following normal distribution
    void gasdev_seq(int &idum, std::vector<double> &vec, int size, double means, double var);

    // generate a sequence following normal distribution with zero average
    void gasdev_seq(int &idum, std::vector<double> &vec, int size, double var);

    // generate a gamma random variable when alpha is larger than 1.0
    double cheng_gamdev(int &idum, const double alpha);

    // generate a random variable follow chi-square distribution
    double chidev(int &idum, const double df);

    // generate a integer between a and b
    int RandAbs(int a, int b, int &seed);

    // Function for get a Chi value of right tail when given a df and prob.
    double chi_val(double df, double prob);

    // Function for get a t value of right tail when given a df and prob.
    double t_val(double df, double prob);

    // Function for get a F value of right tail when given a df and prob.
    double F_val(double df_1, double df_2, double prob);

    // Control the experimental-wise type I error by FDR method
    double ControlFDR(const std::vector<double> &P_Value, double alpha, bool Restrict);
    double ControlFDR_Zou(const std::vector<double> &GenePValue, double FDR);
    double ControlFDR_Storey(std::vector<double> &P_Value, std::vector<double> &Q_Value, double CrtQ, double &FDR);
    double CalcuPi0(std::vector<double> &P_Value, std::vector<double> &Lambda);
    void spline(std::vector<double> &x, std::vector<double> &y, const double yp1, const double ypn, std::vector<double> &y2);
    void splint(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, const double x, double &y);
    std::vector<double> ControlFDR_BH(const std::vector<double> p_value);

    // normal distribution
    double erf(double x);
    double pnorm(double x);
    double dnorm(double x);
    double qnorm_sub(double x, double y);
    double qnorm(double p, bool upper = true);

    // chisq distribution
    double pchisq(double x, double df);
    double qchisq(double q, double df);
    
    // sum of chisq distribution
    double pchisqsum(double x, Eigen::VectorXd lambda);
    double psadd(double x, Eigen::VectorXd lambda);
    double psatt(double x, Eigen::VectorXd lambda);
    double K(double zeta, Eigen::VectorXd &lambda);
    double Kp(double zeta, Eigen::VectorXd &lambda);
    double Kpp(double zeta, Eigen::VectorXd &lambda);
    double Kp_min_x(double zeta, Eigen::VectorXd &lambda, double x);
    double Brents_Kp_min_x(Eigen::VectorXd &lambda, double x, double lowerLimit, double upperLimit, double errorTol);
    std::vector<size_t> sort_re_index(const std::vector<double> &x);
}

#endif
