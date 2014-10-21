#include <Rcpp.h>
#include "MatTypes.h"
#include "SVDsSym.h"

using Rcpp::as;


SEXP do_svds_gen(MatOp *op, SEXP m_scalar_r, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r);


// Main function to calculate truncated SVD for
// dense, real, symmetric matrices
RcppExport SEXP svds_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
BEGIN_RCPP
   
    Rcpp::List params_svds(params_list_r);
    
    int n = as<int>(n_scalar_r);
    int k = as<int>(k_scalar_r);
    int nu = as<int>(nu_scalar_r);
    int nv = as<int>(nv_scalar_r);
    int ncv = as<int>(params_svds["ncv"]);
    double tol = as<double>(params_svds["tol"]);
    int maxitr = as<int>(params_svds["maxitr"]);
    bool rvec = (nu > 0) || (nv > 0);

    MatOp *op = newMatOp(A_mat_r, as<int>(mattype_scalar_r), n, n);
    
    SVDsSym svd(n, k, nu, nv, ncv, op, tol, maxitr);
    svd.compute(rvec);
    SEXP res = svd.extract();

    delete op;

    return res;

END_RCPP
}

// Main function to calculate truncated SVD for
// real general matrices
RcppExport SEXP svds_gen(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                         SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
BEGIN_RCPP

    int m = as<int>(m_scalar_r);
    int n = as<int>(n_scalar_r);
    MatOp *op = newMatOp(A_mat_r, as<int>(mattype_scalar_r), m, n);

    SEXP res = do_svds_gen(op, m_scalar_r, n_scalar_r, k_scalar_r,
                           nu_scalar_r, nv_scalar_r,
                           params_list_r);

    delete op;

    return res;

END_RCPP
}
