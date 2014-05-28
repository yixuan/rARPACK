#ifndef MATOPDGCMATRIX_H
#define MATOPDGCMATRIX_H

#include <RcppEigen.h>
#include "MatOp.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::SparseLU;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix< std::complex<double> > SpCMat;
typedef Eigen::MappedSparseMatrix<double> MapSpMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;

// Operations on "dgCMatrix" class, defined in Matrix package
class MatOp_dgCMatrix : public MatOp
{
private:
    // Sparse matrix structure
    MapSpMat A;
    // Matrix inverse solver
    SparseLU<SpMat> solver;
    SparseLU<SpCMat> csolver;
    VectorXcd cx_vec;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;
public:
    // Constructor
    MatOp_dgCMatrix(SEXP mat_, double sigmar_ = 0, double sigmai_ = 0,
                    bool needSolve_ = false);
    // y_out = A * x_in
    void prod(double *x_in, double *y_out);
    // y_out = A' * x_in
    void tprod(double *x_in, double *y_out);
    // y_out = inv(A - sigma * I) * x_in
    void shiftSolve(double *x_in, double *y_out);
    // Destructor
    virtual ~MatOp_dgCMatrix() {}
};


#endif // MATOPDGCMATRIX_H
