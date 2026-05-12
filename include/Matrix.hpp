#ifndef GCTA2_MATRIX_H
#define GCTA2_MATRIX_H

#include "cpu.h"
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <Logger.h>

static_assert(std::numeric_limits<double>::is_iec559, "Not a supported compiler");

// two step to inverse the matrix

enum INVmethod{
    INV_LLT = 1,
    INV_LU = 2,
    INV_QR = 3,
    INV_FQR = 4,
    INV_SVD = 5,
    INV_ERR = 100
};

template<typename MatrixType>
bool SquareMatrixInverse(MatrixType &A, double &logdet, int &rank, INVmethod &method){
    int n = A.rows();
    bool ret = false;
    switch(method){
        case INV_LLT:
            {
                // Cholesky factorization + inversion via LAPACK dpotrf/dpotri.
                // Both operate in-place on the lower triangle; A is restored on failure
                // so the LU fallback receives the original matrix.
                const MatrixType A_orig = A;
                gcta_blas_int blas_n = static_cast<gcta_blas_int>(n);
                if (gcta_dpotrf(blas_n, A.data(), blas_n) == 0) {
                    // log|det(A)| = 2 * sum(log(diag(L)))  [= sum(log(diag(L)²))]
                    logdet = A.diagonal().array().square().log().sum();
                    if (gcta_dpotri(blas_n, A.data(), blas_n) == 0) {
                        A.template triangularView<Eigen::Upper>() = A.transpose();
                        method = INV_LLT;
                        ret = true;
                        break;
                    }
                }
                A = A_orig;
            }
            [[fallthrough]];
        case INV_LU:
            {
                // LU factorization + inversion via LAPACK dgetrf/dgetri.
                // Uses blocked BLAS3 algorithms — faster than Eigen's PartialPivLU for large n.
                // A is restored on failure so the QR fallback receives the original matrix.
                const MatrixType A_orig = A;
                gcta_blas_int blas_n = static_cast<gcta_blas_int>(n);
                std::vector<gcta_blas_int> ipiv(blas_n);
                if (gcta_dgetrf(blas_n, A.data(), blas_n, ipiv.data()) == 0) {
                    // log|det(A)| = sum(log|diag(U)|); L has unit diagonal
                    logdet = A.diagonal().array().abs().log().sum();
                    if (gcta_dgetri(blas_n, A.data(), blas_n, ipiv.data()) == 0) {
                        method = INV_LU;
                        ret = true;
                        break;
                    }
                }
                A = A_orig;
            }
            [[fallthrough]];
        case INV_QR:
            {
                Eigen::HouseholderQR<MatrixType> qr(A);
                // Compute logdet once instead of absDeterminant() + logAbsDeterminant()
                logdet = qr.logAbsDeterminant();
                if(logdet > std::log(1e-16)){
                    A = qr.solve(MatrixType::Identity(n, n));
                    method = INV_QR;
                    ret = true;
                    break;
                }
            }
            [[fallthrough]];
        case INV_FQR:
            // not necessary 
            // Eigen::HouseholderQR<MatrixType> qr(A);
            {
                Eigen::ColPivHouseholderQR<MatrixType> qr(A);
                if(qr.isInvertible()){
                    logdet = qr.logAbsDeterminant();
                    A = qr.inverse();
                    method = INV_QR;
                    ret = true;
                    break;
                }else{
                    rank = qr.rank();
                    // it will be extreme slow or not accurate
                    if(n > 50000 || 1.0 * rank / n < 0.99){
                        ret = false;
                        break;
                    }
                }
            }
            [[fallthrough]];
        case INV_SVD:
            //Eigen::BDCSVD<MatrixType> svd(A, Eigen::ComputeThinU|Eigen::ComputeThinV);
            ;
        default:
            rank = 0;
            ret = false;
            method = INV_ERR;
    }
    return ret;
}

template <typename T> 
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif //header

