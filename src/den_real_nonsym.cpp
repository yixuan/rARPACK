#include "do_eigs.h"

// Dense matrix-vector product
void den_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    char trans = 'N';
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dgemv)(&trans, &n, &n,
            &alpha, REAL(mat), &n,
            x_in, &one, &zero,
            y_out, &one);
}


RcppExport SEXP den_real_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r,
        SEXP sigmar_scalar_r, SEXP sigmai_scalar_r)
{
BEGIN_RCPP

    return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                   which_string_r, ncv_scalar_r,
                   tol_scalar_r, maxitr_scalar_r,
                   retvec_logical_r,
                   sigmar_scalar_r, sigmai_scalar_r,
                   den_mat_v_prod, NULL);

END_RCPP
}

