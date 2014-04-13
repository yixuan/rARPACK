#include <RcppEigen.h>
#include "ARPACK.h"
#include "EigsSym.h"
#include "MatOpMatrix.h"
#include "MatOpSymMatrix.h"
#include "MatOpDGCMatrix.h"

SEXP do_svds_sym(MatOp *op, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r);

SEXP do_svds_gen(MatOp *op, SEXP m_scalar_r, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r);


// Main function to calculate truncated SVD for
// dense, real, symmetric matrices
RcppExport SEXP den_real_sym_svd(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                 SEXP nu_scalar_r, SEXP nv_scalar_r,
                                 SEXP params_list_r,
                                 SEXP lower_logical_r)
{
BEGIN_RCPP

    char uplo = LOGICAL(lower_logical_r)[0] ? 'L' : 'U';
    MatOpSymMatrix op(A_mat_r, uplo, 0, false);

    return do_svds_sym(&op, n_scalar_r, k_scalar_r,
                       nu_scalar_r, nv_scalar_r,
                       params_list_r);

END_RCPP
}

// Main function to calculate truncated SVD for
// dense, real, general matrices
RcppExport SEXP den_real_gen_svd(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                                 SEXP params_list_r)
{
BEGIN_RCPP

    MatOpMatrix op(A_mat_r, 0, 0, false);

    return do_svds_gen(&op, m_scalar_r, n_scalar_r, k_scalar_r,
                       nu_scalar_r, nv_scalar_r,
                       params_list_r);

END_RCPP
}

// Main function to calculate truncated SVD for
// sparse, real, general matrices
RcppExport SEXP sparse_real_gen_svd(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                                    SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                                    SEXP params_list_r)
{
BEGIN_RCPP

    MatOpDGCMatrix op(A_mat_r, 0, 0, false);

    return do_svds_gen(&op, m_scalar_r, n_scalar_r, k_scalar_r,
                       nu_scalar_r, nv_scalar_r,
                       params_list_r);

END_RCPP
}

