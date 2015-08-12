#include <RcppEigen.h>
#include <SymEigsSolver.h>
#include "matops.h"

using Rcpp::as;

RcppExport SEXP svds_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_svds(params_list_r);

    int n        = as<int>(n_scalar_r);
    int k        = as<int>(k_scalar_r);
    int nu       = as<int>(nu_scalar_r);
    int nv       = as<int>(nv_scalar_r);
    int ncv      = as<int>(params_svds["ncv"]);
    double tol   = as<double>(params_svds["tol"]);
    int maxitr   = as<int>(params_svds["maxitr"]);
    int mattype  = as<int>(mattype_scalar_r);

    MatProd *op = get_mat_prod_op(A_mat_r, n, n, params_list_r, mattype);

    // Prepare initial residuals
    double *init_resid;
    #include "rands.h"
    if(n <= rands_len)
    {
        init_resid = rands;
    } else {
        init_resid = new double[n];
        double *coef_pntr = init_resid;
        for(int i = 0; i < n / rands_len; i++, coef_pntr += rands_len)
        {
            std::copy(rands, rands + rands_len, coef_pntr);
        }
        std::copy(rands, rands + n % rands_len, coef_pntr);
    }

    SymEigsSolver<double, LARGEST_MAGN, MatProd> eigs(op, k, ncv);
    eigs.init(init_resid);
    int nconv = eigs.compute(maxitr, tol);
    if(nconv < k)
        Rcpp::warning("only %d singular values converged, less than k = %d", nconv, k);

    nu = std::min(nu, nconv);
    nv = std::min(nv, nconv);

    Eigen::VectorXd evals = eigs.eigenvalues();
    Eigen::MatrixXd evecs = eigs.eigenvectors();

    Rcpp::NumericVector d(nconv);
    Rcpp::NumericMatrix u(n, nu), v(n, nv);

    // Copy evals to d
    std::copy(evals.data(), evals.data() + nconv, d.begin());

    // Copy evecs to u and v with the specified number of columns
    std::copy(evecs.data(), evecs.data() + nu * n, u.begin());
    std::copy(evecs.data(), evecs.data() + nv * n, v.begin());

    // We need to make sure that singular values are nonnegative,
    // so move the sign to v.
    for(int i = 0; i < nconv; i++)
    {
        if(d[i] < 0)
        {
            d[i] = -d[i];
            if(i < nv)
            {
                double *ptr = &v(0, i);
                std::transform(ptr, ptr + n, ptr, std::negate<double>());
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("d")     = d,
        Rcpp::Named("u")     = u,
        Rcpp::Named("v")     = v,
        Rcpp::Named("niter") = eigs.num_iterations(),
        Rcpp::Named("nops")  = eigs.num_operations() * 2
    );

    END_RCPP
}
