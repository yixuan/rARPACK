#include "MatOp_dgCMatrix.h"

MatOp_dgCMatrix::MatOp_dgCMatrix(SEXP mat_, double sigmar_, double sigmai_,
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

    if(m != n) return;

    // If sigma is real
    if(fabs(sigmai) < 1e-17)
    {
        // Create a sparse idendity matrix
        SpMat I(n, n);
        I.setIdentity();

        // Sparse LU decomposition
        solver.compute(A - sigmar * I);

    } else {
        SpCMat cA = A.cast< std::complex<double> >();
        
        // Create a sparse identity matrix (1 + 0i on diagonal)
        SpCMat I(n, n);
        I.setIdentity();
        
        // Sparse LU decomposition
        csolver.compute(cA - std::complex<double>(sigmar, sigmai) * I);

        cx_vec.resize(n);
        cx_vec.setZero();
    }
}

void MatOp_dgCMatrix::prod(double *x_in, double *y_out)
{
    new (&x_vec) MapVec(x_in, n);
    new (&y_vec) MapVec(y_out, m);
    y_vec = A * x_vec;
}

void MatOp_dgCMatrix::tprod(double *x_in, double *y_out)
{
    new (&x_vec) MapVec(x_in, m);
    new (&y_vec) MapVec(y_out, n);
    y_vec = A.transpose() * x_vec;
}

void MatOp_dgCMatrix::shiftSolve(double *x_in, double *y_out)
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
