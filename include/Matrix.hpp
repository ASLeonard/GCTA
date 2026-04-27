#ifndef GCTA2_MATRIX_H
#define GCTA2_MATRIX_H

#include "cpu.h"
#include <Eigen/Eigen>
#include <iostream>
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

/*
template<typename MatrixType>
bool _LLT(MatrixType &A, double &logdet){
    //MatrixType::Scalar * vi = A.data();
    auto * vi = A.data();
    Eigen::Matrix<typename MatrixType::Scalar, Eigen::Dynamic, 1> diag = A.diagonal();
    //auto diag = A.diagonal();

    gcta_blas_int info = 0;
    gcta_blas_int cols = (gcta_blas_int)A.cols();
    char uplo = 'L';
    dpotrf(&uplo, &cols, vi, &cols, &info);
    if(info == 0){
        logdet = A.diagonal().array().square().log().sum();
        dpotri(&uplo, &cols, vi, &cols, &info);
        if(info == 0){
            A.template triangularView<Eigen::Upper>() = A.transpose();
            return true;
        }
    }

    A.template triangularView<Eigen::Lower>() = A.transpose();
    A.diagonal() = diag;
    return false;
}


// LU factorization + inversion via LAPACK dgetrf/dgetri.
// Uses blocked BLAS3 algorithms — faster than Eigen's PartialPivLU for large n.
// Restores A to its original state on failure so downstream fallbacks see the
// unmodified matrix (dgetrf overwrites A in-place).
template<typename MatrixType>
bool _LU(MatrixType &A, double &logdet) {
    MatrixType A_orig = A;
    auto *vi = A.data();
    gcta_blas_int n = (gcta_blas_int)A.cols();
    gcta_blas_int info = 0;
    std::vector<gcta_blas_int> ipiv(n);

    dgetrf(&n, &n, vi, &n, ipiv.data(), &info);
    if (info != 0) { A = A_orig; return false; }

    // log|det(A)| = sum(log|diag(U)|); L has unit diagonal so doesn't contribute
    logdet = A.diagonal().array().abs().log().sum();

    // Workspace query, then in-place inversion
    gcta_blas_int lwork = -1;
    double work_query = 0.0;
    dgetri(&n, vi, &n, ipiv.data(), &work_query, &lwork, &info);
    lwork = (gcta_blas_int)work_query;
    std::vector<double> work(lwork);
    dgetri(&n, vi, &n, ipiv.data(), work.data(), &lwork, &info);
    if (info != 0) { A = A_orig; return false; }
    return true;
}
*/

template<typename MatrixType>
bool _LLT(MatrixType &A, double &logdet) {
    auto diag = A.diagonal().eval();

    Eigen::LLT<MatrixType> llt(A);
    if (llt.info() != Eigen::Success) {
        A.diagonal() = diag;
        return false;
    }

    logdet = llt.matrixL().nestedExpression().diagonal().array().log().sum() * 2.0;
    A = llt.solve(MatrixType::Identity(A.rows(), A.cols()));
    return true;
}

template<typename MatrixType>
bool _LU(MatrixType &A, double &logdet) {
    Eigen::PartialPivLU<MatrixType> lu(A);

    logdet = lu.matrixLU().diagonal().array().abs().log().sum();
    if (!std::isfinite(logdet)) {
        return false;
    }

    A = lu.inverse();
    return true;
}

template<typename MatrixType>
bool SquareMatrixInverse(MatrixType &A, double &logdet, int &rank, INVmethod &method){
    int n = A.rows();
    bool ret = false;
    switch(method){
        case INV_LLT:
            {
                if(_LLT(A, logdet)){
                    method = INV_LLT;
                    ret = true;
                    break;
                }
            }
            [[fallthrough]];
        case INV_LU:
            {
                // LAPACK dgetrf/dgetri: blocked BLAS3, faster than Eigen PartialPivLU.
                // _LU restores A on failure so INV_QR fallback gets the original matrix.
                if(_LU(A, logdet)){
                    method = INV_LU;
                    ret = true;
                    break;
                }
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

