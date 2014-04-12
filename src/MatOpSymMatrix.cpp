#include "MatOpSymMatrix.h"

const double MatOpSymMatrix::BLAS_alpha = 1.0;
const int MatOpSymMatrix::BLAS_one = 1;
const double MatOpSymMatrix::BLAS_zero = 0.0;

MatOpSymMatrix::MatOpSymMatrix(SEXP mat_, char uplo_, double sigma_,
                               bool needSolve_) :
    A_pntr(REAL(mat_)), A(as<MapMat>(mat_)), uplo(uplo_),
    x_vec(NULL, n), y_vec(NULL, n)
{
    m = n = A.rows();
    sigmar = sigma_;
    sigmai = 0;
    canTprod = true;
    canSolve = needSolve_;

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

void MatOpSymMatrix::prod(double *x_in, double *y_out)
{
    F77_CALL(dsymv)(&uplo, &n,
            &BLAS_alpha, A_pntr, &n,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void MatOpSymMatrix::tprod(double *x_in, double *y_out)
{
    F77_CALL(dsymv)(&uplo, &n,
            &BLAS_alpha, A_pntr, &n,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void MatOpSymMatrix::shiftSolve(double *x_in, double *y_out)
{
    if(!canSolve)
        Rcpp::stop("this matrix doesn't support solving linear equation");

    new (&x_vec) MapVec(x_in, n);
    new (&y_vec) MapVec(y_out, n);
    y_vec = solver.solve(x_vec);
}
