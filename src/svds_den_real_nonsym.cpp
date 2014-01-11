#include <RcppEigen.h>
#include "do_eigs.h"

// Calculate y_out = Mat * x_in
static void den_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                           int m, int n, void *data)
{
    char trans = 'N';
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dgemv)(&trans, &m, &n,
            &alpha, REAL(mat), &m,
            x_in, &one, &zero,
            y_out, &one);
}

// Calculate y_out = Mat' * x_in
static void den_mat_t_v_prod(SEXP mat, double *x_in, double *y_out,
                             int m, int n, void *data)
{
    char trans = 'T';
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dgemv)(&trans, &m, &n,
            &alpha, REAL(mat), &m,
            x_in, &one, &zero,
            y_out, &one);
}


// Main function to calculate truncated SVD for
// dense, real, nonsymmetric matrices
RcppExport SEXP den_real_nonsym_svd(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                                    SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                                    SEXP params_list_r)
{
BEGIN_RCPP

    return do_svds_nonsym(A_mat_r, m_scalar_r, n_scalar_r,
                          k_scalar_r, nu_scalar_r, nv_scalar_r,
                          params_list_r, den_mat_v_prod, den_mat_t_v_prod,
                          NULL);

END_RCPP
}

