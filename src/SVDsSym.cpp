#include <RcppEigen.h>
#include "do_eigs.h"

using Rcpp::as;
using Rcpp::wrap;

// For symmetric matrices, svds() is equivalent to eigs(),
// so many steps can be simplified.
//
// Main function to calculate real, symmetric SVD
SEXP do_svds_sym(SEXP A_mat_r, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r, Mvfun mat_v_prod,
                 void *data)
{
BEGIN_RCPP

    // We will use do_eigs_sym() to calculate the results,
    // but the parameters list is slightly different between
    // eigs() and svds().
    //
    // We create a new list params_eigs to be passed to
    // do_eigs_sym().

    // Retrieve parameters
    Rcpp::List params_svds(params_list_r);
    int n = INTEGER(n_scalar_r)[0];
    int k = INTEGER(k_scalar_r)[0];
    // Whether to calculate singular vectors or not.
    int nu = INTEGER(nu_scalar_r)[0];
    int nv = INTEGER(nv_scalar_r)[0];
    bool rvec = (nu > 0) | (nv > 0);
    // Create list for do_eigs_sym()
    Rcpp::List params_eigs =
        Rcpp::List::create(Rcpp::Named("which") = wrap("LM"),
                           Rcpp::Named("ncv") = params_svds["ncv"],
                           Rcpp::Named("tol") = params_svds["tol"],
                           Rcpp::Named("maxitr") = params_svds["maxitr"],
                           Rcpp::Named("retvec") = wrap(rvec),
                           Rcpp::Named("sigma") = wrap(0.0),
                           Rcpp::Named("workmode") = wrap(1L));
    // Result from do_eigs_sum()
    SEXP eigs_res = do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                                params_eigs, mat_v_prod, data);

    // Result list
    Rcpp::List ret(eigs_res);

    int nconv = as<int>(ret["nconv"]);
    if (nconv < k)
        ::Rf_warning("only %d singular values converged, less than k", nconv);
    nu = nu > nconv ? nconv : nu;
    nv = nv > nconv ? nconv : nv;

    // Currently ret has components of values, vectors, nconv
    // and niter, while what we need is d, u, v, nconv and niter,
    // so we first insert the v matrix into the list, and then
    // change the list names.
    if (!rvec)
    {
        ret.insert(2, R_NilValue);
    } else {
        // At least one of nu and nv are not zero
        if (nu != 0)
        {
            Rcpp::NumericMatrix u = ret["vectors"];
            if (nv == 0)
                // nu != 0, nv == 0
                ret.insert(2, R_NilValue);
            else {
                // nu != 0, nv != 0
                Rcpp::NumericMatrix v(n, nv);
                // In SVD, V == U
                std::copy(u.begin(),
                          u.begin() + nv * n,
                          v.begin());
                ret.insert(2, v);
            }
            // u.erase(start, end) removes u[start <= i < end]
            u.erase(nu * n, nconv * n);
        } else {
            // nu == 0, nv != 0
            Rcpp::NumericMatrix u = ret["vectors"];
            Rcpp::NumericMatrix v(n, nv);
            std::copy(u.begin(),
                      u.begin() + nv * n,
                      v.begin());
            ret["vectors"] = R_NilValue;
            ret.insert(2, v);
        }
    }
    ret.names() = Rcpp::CharacterVector::create("d", "u", "v",
                                                "nconv", "niter");

    return ret;

END_RCPP
}

