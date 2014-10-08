#ifndef MATOP_FUNCTION_H
#define MATOP_FUNCTION_H

#include <RcppEigen.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL
#include "MatOp.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::PartialPivLU;

typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;

// Operations on "matrix" class, defined in base R
class MatOp_function : public MatOp
{
private:
    // Pointer to the A matrix
    Rcpp::Function FUN_function_;
    // Some constants for BLAS
    static const char BLAS_notrans;
    static const char BLAS_trans;
    static const double BLAS_alpha;
    static const int BLAS_one;
    static const double BLAS_zero;
    // Matrix inverse solver
    // PartialPivLU<MatrixXd> solver;
    // PartialPivLU<MatrixXcd> csolver;
    // VectorXcd cx_vec;
    // Mapped vector
    // MapVec x_vec;
    // MapVec y_vec;
    // Function called by constructors
    void init();
public:
    // Constructor
    MatOp_function(Rcpp::Function FUN_function_r, int m_, int n_,
                 bool needSolve_ = false);
    // y_out = FUN(x_in) 
    void prod(double *x_in, double *y_out);
    // y_out = FUN(x_in) suppose it's a symetric map at first
    void tprod(double *x_in, double *y_out);
    // // y_out = inv(A - sigma * I) * x_in
    // void shiftSolve(double *x_in, double *y_out);
    // Destructor
    virtual ~MatOp_function() {}
};


#endif // MATOP_MATRIX_H
