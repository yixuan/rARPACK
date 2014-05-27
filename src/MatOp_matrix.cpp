#include "MatOp_matrix.h"

const char MatOp_matrix::BLAS_notrans = 'N';
const char MatOp_matrix::BLAS_trans = 'T';
const double MatOp_matrix::BLAS_alpha = 1.0;
const int MatOp_matrix::BLAS_one = 1;
const double MatOp_matrix::BLAS_zero = 0.0;

MatOp_matrix::MatOp_matrix(SEXP mat_, double sigmar_, double sigmai_,
                          bool needSolve_) :
    A_pntr(REAL(mat_)), A(as<MapMat>(mat_)),
    x_vec(NULL, n), y_vec(NULL, n)
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
        // Subtract the diagonal elements by sigma, i.e., A - sigma * I
        for(int i = 0; i < n; i++)
        {
            A(i, i) -= sigmar;
        }
        // LU decomposition
        solver.compute(A);
        // Change A back
        for(int i = 0; i < n; i++)
        {
            A(i, i) += sigmar;
        }
    } else {
        MatrixXcd cA = A.cast< std::complex<double> >();
        // Subtract the diagonal elements by sigma, i.e., A - sigma * I
        for(int i = 0; i < n; i++)
        {
            cA(i, i) -= std::complex<double>(sigmar, sigmai);
        }
        // LU decomposition
        csolver.compute(cA);

        cx_vec.resize(n);
        cx_vec.setZero();
    }
}

void MatOp_matrix::prod(double *x_in, double *y_out)
{
    F77_CALL(dgemv)(&BLAS_notrans, &m, &n,
            &BLAS_alpha, A_pntr, &m,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void MatOp_matrix::tprod(double *x_in, double *y_out)
{
    F77_CALL(dgemv)(&BLAS_trans, &m, &n,
            &BLAS_alpha, A_pntr, &m,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void MatOp_matrix::shiftSolve(double *x_in, double *y_out)
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
