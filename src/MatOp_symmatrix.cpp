#include "MatOp_symmatrix.h"

const double MatOp_symmatrix::BLAS_alpha = 1.0;
const int MatOp_symmatrix::BLAS_one = 1;
const double MatOp_symmatrix::BLAS_zero = 0.0;

MatOp_symmatrix::MatOp_symmatrix(SEXP mat_, int n_, char uplo_, double sigma_,
                                 bool needSolve_) :
    MatOp(n_, n_, sigma_, 0.0, true, needSolve_),
    A_pntr(REAL(mat_)), A(A_pntr, n, n), uplo(uplo_),
    x_vec(NULL, n), y_vec(NULL, n)
{
    if(!needSolve_)  return;

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
}

void MatOp_symmatrix::prod(double *x_in, double *y_out)
{
    F77_CALL(dsymv)(&uplo, &n,
            &BLAS_alpha, A_pntr, &n,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void MatOp_symmatrix::tprod(double *x_in, double *y_out)
{
    F77_CALL(dsymv)(&uplo, &n,
            &BLAS_alpha, A_pntr, &n,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void MatOp_symmatrix::shiftSolve(double *x_in, double *y_out)
{
    if(!canSolve)
        Rcpp::stop("this matrix doesn't support solving linear equation");

    new (&x_vec) MapVec(x_in, n);
    new (&y_vec) MapVec(y_out, n);
    y_vec = solver.solve(x_vec);
}
