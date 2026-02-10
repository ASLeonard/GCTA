/*
 * Implementations of the extended function for EIGEN lib
 *
 * 2014 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "eigen_func.h"

void eigen_func::rank(Eigen::VectorXf &x, Eigen::VectorXf &rank)
{
    int size = x.size();
    if(size < 1) return;

    int i = 0;
    rank.resize(size);
    std::vector<float> x_sort(size);
    for(i = 0; i < size; i++) x_sort[i] = x[i];
    std::stable_sort(x_sort.begin(), x_sort.end());
    std::map<double, int> value_indx;
    for(i = 0; i < size; i++) value_indx.insert(std::pair<double,int>(x_sort[i], i+1));
    std::map<double, int>::iterator iter;
    for(i = 0; i < size; i++){
        iter = value_indx.find(x[i]);
        rank[i] = iter->second;
    }
}

void eigen_func::inverse_norm_rank_transform(Eigen::VectorXf &x)
{
    int size = x.size();
    if(size < 1) return;

    Eigen::VectorXf rank;
    eigen_func::rank(x, rank);

    int i = 0;
    float size_f = (float) size;
    for(i = 0; i < size; i++) x[i] = StatFunc::qnorm((rank[i] - 0.5) / size_f);
}
