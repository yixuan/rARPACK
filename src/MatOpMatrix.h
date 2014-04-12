#ifndef MATOPMATRIX_H
#define MATOPMATRIX_H

#include <RcppEigen.h>
#include "MatOp.h"
#include "do_eigs.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::PartialPivLU;

typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;

class MatOpMatrix : public MatOp
{
private:
    // Pointer to the A matrix
    double *A_pntr;
    // Map mat_ to Eigen matrix
    MapMat A;
    // Some constants for BLAS
    static const char BLAS_notrans;
    static const char BLAS_trans;
    static const double BLAS_alpha;
    static const int BLAS_one;
    static const double BLAS_zero;
    // Matrix inverse solver
    PartialPivLU<MatrixXd> solver;
    PartialPivLU<MatrixXcd> csolver;
    VectorXcd cx_vec;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;
public:
    // Constructor
    MatOpMatrix(SEXP mat_, double sigmar_ = 0, double sigmai_ = 0,
                bool needSolve_ = false);
    // y_out = A * x_in
    void prod(double *x_in, double *y_out);
    // y_out = A' * x_in
    void tprod(double *x_in, double *y_out);
    // y_out = inv(A - sigma * I) * x_in
    void shiftSolve(double *x_in, double *y_out);
    // Destructor
    ~MatOpMatrix() {}
};


#endif // MATOPMATRIX_H
