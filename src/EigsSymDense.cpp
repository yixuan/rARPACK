#include "EigsSymDense.h"

const double EigsSymDense::BLAS_alpha = 1.0;
const int EigsSymDense::BLAS_one = 1;
const double EigsSymDense::BLAS_zero = 0.0;

EigsSymDense::EigsSymDense(int n_, int nev_, int ncv_,
                           const string & which_, int workmode_,
                           double sigmar_,
                           char bmat_, double tol_, int maxitr_) :
    EigsSym(n_, nev_, ncv_, which_, workmode_,
            sigmar_, bmat_, tol_, maxitr_),
    x_vec(NULL, n, n),
    y_vec(NULL, n, n)
{
    A_pntr = NULL;
    uplo = 'L';
}

EigsSymDense::~EigsSymDense()
{
    delete [] workl;
    delete [] workd;
    delete [] resid;
}

void EigsSymDense::BindMatrix(SEXP mat_, char uplo_)
{
    A_pntr = REAL(mat_);
    uplo = uplo_;

    // If we don't need shift-and-invert mode
    if (workmode == 1)
    {
        MatrixLinked(true);
        return;
    }

    // Map mat_ to Eigen matrix
    MapMat A(A_pntr, n, n);
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
    MatrixLinked(true);
}

void EigsSymDense::BindMatrix(SEXP mat_)
{
    BindMatrix(mat_, 'L');
}

void EigsSymDense::MultVector(double *x_in, double *y_out)
{
    F77_CALL(dsymv)(&uplo, &n,
            &BLAS_alpha, A_pntr, &n,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void EigsSymDense::MultVectorShift(double *x_in, double *y_out)
{
    new (&x_vec) MapVec(x_in, n);
    new (&y_vec) MapVec(y_out, n);
    y_vec = solver.solve(x_vec);
}

