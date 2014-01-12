#include <RcppEigen.h>
#include "do_eigs.h"

// Dense matrix-vector product
// This function uses BLAS to calculte y_out = Mat * x_in
// where Mat is a symmetric matrix
static void den_sym_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                               int m, int n, void *data)
{
    // *data is eigher 'L' or 'U'
    char *uplo = (char *) data;
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dsymv)(uplo, &n,
             &alpha, REAL(mat), &n,
             x_in, &one, &zero,
             y_out, &one);
}

// Main function to calculate truncated SVD for
// dense, real, symmetric matrices
RcppExport SEXP den_real_sym_svd(SEXP A_mat_r, SEXP n_scalar_r,
                                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                                 SEXP params_list_r,
                                 SEXP lower_logical_r)
{
BEGIN_RCPP

    char uplo = LOGICAL(lower_logical_r)[0] ? 'L' : 'U';

    return do_svds_sym(A_mat_r, n_scalar_r,
                       k_scalar_r, nu_scalar_r, nv_scalar_r,
                       params_list_r, den_sym_mat_v_prod,
                       &uplo);

END_RCPP
}

