#ifndef MATPROD_MATRIX_H
#define MATPROD_MATRIX_H

#include <Rcpp.h>
#include <R_ext/BLAS.h>  // for BLAS and F77_CALL
#include "MatProd.h"

class MatProd_matrix: public MatProd
{
private:
    const double* mat_ptr;
    const int nrow;
    const int ncol;
    const double BLAS_alpha;
    const int BLAS_one;
    const double BLAS_zero;

public:
    MatProd_matrix(SEXP mat_, const int nrow_, const int ncol_) :
        mat_ptr(REAL(mat_)),
        nrow(nrow_),
        ncol(ncol_),
        BLAS_alpha(1.0),
        BLAS_one(1),
        BLAS_zero(0.0)
    {}

    int rows() { return nrow; }
    int cols() { return ncol; }

    // y_out = A * x_in
    void perform_op(double *x_in, double *y_out)
    {
        F77_CALL(dgemv)("N", &nrow, &ncol,
                        &BLAS_alpha, mat_ptr, &nrow,
                        x_in, &BLAS_one, &BLAS_zero,
                        y_out, &BLAS_one);
    }
};


#endif // MATPROD_MATRIX_H
