/*
 * Implementations of the commonly-used functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "CommFunc.h"
#include "Logger.h"
#include <cmath>
#include <numeric>

double CommFunc::median(std::vector<double> x)
{
    size_t n = x.size();
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    if (n == 1) return x[0];
    size_t mid = n / 2;
    std::nth_element(x.begin(), x.begin() + mid, x.end());
    if (n % 2 == 1) return x[mid];
    double upper = x[mid];
    double lower = *std::max_element(x.begin(), x.begin() + mid);
    return (lower + upper) / 2.0;
}

double CommFunc::quantile(const double* data, size_t n, double prob) {
    if (prob < 0 || prob > 1) LOGGER.e(0, "requested quantile probability is invalid");
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    double index = prob * (n - 1);
    size_t below = std::floor(index), above = std::ceil(index);
    if (below == above) return data[above];
    return (above - index) * data[below] + (index - below) * data[above];
}

double CommFunc::quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob) {
    return quantile(vals.data(), vals.size(), prob);
}

double CommFunc::quantile(const std::vector<double> &vals, double prob) {
    return quantile(vals.data(), vals.size(), prob);
}

double CommFunc::var(const std::vector<double> &x)
{
    if (x.size() <= 1) return 0.0;
    double mu = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    double s2 = std::accumulate(x.begin(), x.end(), 0.0,
        [mu](double acc, double v) { return acc + (v - mu) * (v - mu); });
    return s2 / (x.size() - 1);
}

//Note: float tolerance is used here, which is more relaxed than double tolerance,
//because the FDR values are usually close to 0 or 1.
bool CommFunc::FloatEqual(double lhs, double rhs)
{
    return std::abs(lhs - rhs) < std::numeric_limits<float>::epsilon();
}

bool CommFunc::FloatNotEqual(double lhs, double rhs)
{
    return std::abs(lhs - rhs) >= std::numeric_limits<float>::epsilon();
}
