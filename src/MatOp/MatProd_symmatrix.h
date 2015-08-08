#ifndef MATPROD_SYMMATRIX_H
#define MATPROD_SYMMATRIX_H

#include <Rcpp.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL

class MatProd_symmatrix
{
private:
    const double* mat_ptr;
    const int n;
    const char uplo;
    const double BLAS_alpha;
    const int BLAS_one;
    const double BLAS_zero;

public:
    MatProd_symmatrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        mat_ptr(REAL(mat_)),
        n(nrow_),
        uplo(uplo_),
        BLAS_alpha(1.0),
        BLAS_one(1),
        BLAS_zero(0.0)
    {}

    int rows() { return n; }
    int cols() { return n; }

    // y_out = A * x_in
    void perform_op(double *x_in, double *y_out)
    {
        F77_CALL(dsymv)(&uplo, &n,
                        &BLAS_alpha, mat_ptr, &n,
                        x_in, &BLAS_one, &BLAS_zero,
                        y_out, &BLAS_one);
    }
};


#endif // MATPROD_SYMMATRIX_H
