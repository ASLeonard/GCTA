/*
 * Interface to the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#ifndef _COMMFUNC_H
#define _COMMFUNC_H

#include <cstdlib>
#include <limits>
#include <vector>
#include <algorithm>
#include <Eigen/StdVector>


namespace CommFunc
{
    double median(std::vector<double> x);
    double quantile(const double* data, size_t n, double prob);
    double quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob);
    double quantile(const std::vector<double> &vals, double prob);
    double var(const std::vector<double> &x);
    bool FloatEqual(double lhs, double rhs);
    bool FloatNotEqual(double lhs, double rhs);
}

#endif
