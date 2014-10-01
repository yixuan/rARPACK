#ifndef MATOP_SYMMATRIX_H
#define MATOP_SYMMATRIX_H

#include <RcppEigen.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL
#include "MatOp.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::PartialPivLU;

typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;

class MatOp_symmatrix : public MatOp
{
private:
    // Pointer to the A matrix
    double *A_pntr;
    // Map mat_ to Eigen matrix
    MapMat A;
    // Which part of A should be used
    // 'L' or 'U'
    char uplo;
    // Some constants for BLAS
    static const double BLAS_alpha;
    static const int BLAS_one;
    static const double BLAS_zero;
    // Matrix inverse solver
    PartialPivLU<MatrixXd> solver;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;
public:
    // Constructor
    MatOp_symmatrix(SEXP mat_, int n_, char uplo_ = 'L',
                    double sigma_ = 0,
                    bool needSolve_ = false);
    // y_out = A * x_in
    void prod(double *x_in, double *y_out);
    // y_out = A' * x_in
    void tprod(double *x_in, double *y_out);
    // y_out = inv(A - sigma * I) * x_in
    void shiftSolve(double *x_in, double *y_out);
    // Destructor
    virtual ~MatOp_symmatrix() {}
};


#endif // MATOP_SYMMATRIX_H
