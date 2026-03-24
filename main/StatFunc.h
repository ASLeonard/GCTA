/*
 * Interface to the statistical functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 * Modernized 2026 with C++20/23/26 features
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
#include <ranges>
#include <span>
#include <numbers>
#include <iostream>
#include <random>
#include "CommFunc.h"
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

namespace StatFunc
{
    ////////// P-value Calculation Functions Start ////////////////
    // if two_tail=true, return two tail test probability;
    // otherwise, return one tail test probability;
    [[nodiscard]] double t_prob(double df, double t_value, bool two_tail = true);

    [[nodiscard]] double F_prob(double df_1, double df_2, double F_value);

    // log_p=true returns ln(p) instead of p, avoiding underflow for extreme chi-squared values
    [[nodiscard]] double chi_prob(double df, double chi_sqr_val, bool log_p = false);
    ////////// P-value Calculatiion Functions End ////////////////

    ///////// Random Number Generation Functions Start ////////
    // (0, 1) uniform distribution generator — idum is accepted for API compatibility but ignored
    [[nodiscard]] inline double uniform_random_distribution(int & /*idum*/) {
        static std::random_device rd;
        static std::mt19937 rng(rd());
        static std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(rng);
    }

    // (a, b) uniform distribution generator
    [[nodiscard]] double UniformDev(double a, double b, int &idum);

    // normal distribution generator
    [[nodiscard]] double gasdev(int &idum);

    // generate a sequence following normal distribution
    void gasdev_seq(int &idum, std::vector<double> &vec, int size, double means, double var);

    // generate a sequence following normal distribution with zero average
    void gasdev_seq(int &idum, std::vector<double> &vec, int size, double var);

    // generate a gamma random variable when alpha is larger than 1.0
    [[nodiscard]] double cheng_gamdev(int &idum, double alpha);

    // generate a random variable follow chi-square distribution
    [[nodiscard]] double chidev(int &idum, double df);

    // generate an integer between a and b
    [[nodiscard]] int RandAbs(int a, int b, int &seed);

    // Function for get a Chi value of right tail when given a df and prob.
    [[nodiscard]] double chi_val(double df, double prob);

    // Function for get a t value of right tail when given a df and prob.
    [[nodiscard]] double t_val(double df, double prob);

    // Function for get a F value of right tail when given a df and prob.
    [[nodiscard]] double F_val(double df_1, double df_2, double prob);

    // Control the experimental-wise type I error by FDR method
    [[nodiscard]] double ControlFDR(std::span<const double> P_Value, double alpha, bool Restrict);
    [[nodiscard]] double ControlFDR_Zou(std::span<const double> GenePValue, double FDR);
    [[nodiscard]] double ControlFDR_Storey(std::vector<double> &P_Value, std::vector<double> &Q_Value, double CrtQ, double &FDR);
    [[nodiscard]] double CalcuPi0(std::vector<double> &P_Value, std::vector<double> &Lambda);
    void spline(std::vector<double> &x, std::vector<double> &y, double yp1, double ypn, std::vector<double> &y2);
    void splint(std::vector<double> &xa, std::vector<double> &ya, std::vector<double> &y2a, double x, double &y);
    [[nodiscard]] std::vector<double> ControlFDR_BH(std::span<const double> p_value);

    // normal distribution (upper-tail CDF, PDF, and quantile)
    [[nodiscard]] double pnorm(double x);
    [[nodiscard]] double dnorm(double x);
    [[nodiscard]] double qnorm(double p, bool upper = true);

    // chisq distribution
    // log_p=true returns ln(p) instead of p, avoiding underflow for extreme chi-squared values
    [[nodiscard]] double pchisq(double x, double df, bool log_p = false);
    [[nodiscard]] double qchisq(double q, double df);

    // sum of chisq distribution
    [[nodiscard]] double pchisqsum(double x, Eigen::VectorXd lambda);
    [[nodiscard]] double psadd(double x, Eigen::VectorXd lambda);
    [[nodiscard]] double psatt(double x, Eigen::VectorXd lambda);
    [[nodiscard]] double K(double zeta, Eigen::VectorXd &lambda);
    [[nodiscard]] double Kp(double zeta, Eigen::VectorXd &lambda);
    [[nodiscard]] double Kpp(double zeta, Eigen::VectorXd &lambda);
    [[nodiscard]] double Kp_min_x(double zeta, Eigen::VectorXd &lambda, double x);
    [[nodiscard]] double Brents_Kp_min_x(Eigen::VectorXd &lambda, double x, double lowerLimit, double upperLimit, double errorTol);
    [[nodiscard]] std::vector<size_t> sort_re_index(std::span<const double> x);
}

#endif
