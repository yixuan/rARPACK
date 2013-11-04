#include "do_eigs.h"

// Dense matrix-vector product
void den_sym_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    char *uplo = (char *) data;
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dsymv)(uplo, &n,
            &alpha, REAL(mat), &n,
            x_in, &one, &zero,
            y_out, &one);
}


RcppExport SEXP den_real_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r, SEXP sigma_scalar_r,
        SEXP workmode_scalar_r,
        SEXP lower_logical_r)
{
BEGIN_RCPP
    char uplo = LOGICAL(lower_logical_r)[0] ? 'L' : 'U';
    return do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                   which_string_r, ncv_scalar_r,
                   tol_scalar_r, maxitr_scalar_r,
                   retvec_logical_r, sigma_scalar_r,
                   workmode_scalar_r,
                   den_sym_mat_v_prod, &uplo);

END_RCPP
}

