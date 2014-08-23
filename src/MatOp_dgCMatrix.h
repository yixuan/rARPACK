#ifndef MATOPDGCMATRIX_H
#define MATOPDGCMATRIX_H

#include <RcppEigen.h>
#include "MatOp.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::SparseMatrix;
using Eigen::SparseLU;

template<int Storage>
class MatOp_sparseMatrix : public MatOp
{

typedef Eigen::SparseMatrix<double, Storage> SpMat;
typedef Eigen::SparseMatrix< std::complex<double>, Storage > SpCMat;
typedef Eigen::MappedSparseMatrix<double, Storage> MapSpMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;

private:
    // Sparse matrix structure
    MapSpMat A;
    // Matrix inverse solver
    SparseLU< SparseMatrix<double, Eigen::ColMajor> > solver;
    SparseLU< SparseMatrix<std::complex<double>, Eigen::ColMajor> > csolver;
    VectorXcd cx_vec;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;
public:
    // Constructor
    MatOp_sparseMatrix(SEXP mat_, double sigmar_ = 0, double sigmai_ = 0,
                       bool needSolve_ = false);
    // y_out = A * x_in
    void prod(double *x_in, double *y_out);
    // y_out = A' * x_in
    void tprod(double *x_in, double *y_out);
    // y_out = inv(A - sigma * I) * x_in
    void shiftSolve(double *x_in, double *y_out);
    // Destructor
    virtual ~MatOp_sparseMatrix() {}
};


template<int Storage>
MatOp_sparseMatrix<Storage>::MatOp_sparseMatrix(SEXP mat_, double sigmar_, double sigmai_,
                                                bool needSolve_) :
    A(as<MapSpMat>(mat_)),
    x_vec(NULL, 1), y_vec(NULL, 1)
{
    m = A.rows();
    n = A.cols();
    sigmar = sigmar_;
    sigmai = sigmai_;
    canTprod = true;
    canSolve = needSolve_;

    if(!needSolve_)  return;

    if(m != n)  return;

    // If sigma is real
    if(fabs(sigmai) < 1e-17)
    {
        // Create a sparse idendity matrix
        SpMat I(n, n);
        I.setIdentity();

        // Sparse LU decomposition
        solver.compute(A - sigmar * I);

    } else {
        SpCMat cA = A.template cast< std::complex<double> >();
        
        // Create a sparse identity matrix (1 + 0i on diagonal)
        SpCMat I(n, n);
        I.setIdentity();
        
        // Sparse LU decomposition
        csolver.compute(cA - std::complex<double>(sigmar, sigmai) * I);

        cx_vec.resize(n);
        cx_vec.setZero();
    }
}

template<int Storage>
void MatOp_sparseMatrix<Storage>::prod(double *x_in, double *y_out)
{
    new (&x_vec) MapVec(x_in, n);
    new (&y_vec) MapVec(y_out, m);
    y_vec = A * x_vec;
}

template<int Storage>
void MatOp_sparseMatrix<Storage>::tprod(double *x_in, double *y_out)
{
    new (&x_vec) MapVec(x_in, m);
    new (&y_vec) MapVec(y_out, n);
    y_vec = A.transpose() * x_vec;
}

template<int Storage>
void MatOp_sparseMatrix<Storage>::shiftSolve(double *x_in, double *y_out)
{
    if(m != n)
        Rcpp::stop("matrix must be square");
    if(!canSolve)
        Rcpp::stop("this matrix doesn't support solving linear equation");

    if(fabs(sigmai) < 1e-17)
    {
        new (&x_vec) MapVec(x_in, n);
        new (&y_vec) MapVec(y_out, n);
        y_vec = solver.solve(x_vec);
    } else {
        cx_vec.real() = MapVec(x_in, n);
        new (&y_vec) MapVec(y_out, n);
        y_vec = csolver.solve(cx_vec).real();
    }
}

// Operations on "dgCMatrix" class, defined in Matrix package
typedef MatOp_sparseMatrix<Eigen::ColMajor> MatOp_dgCMatrix;


#endif // MATOPDGCMATRIX_H
