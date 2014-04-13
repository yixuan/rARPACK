#include <RcppEigen.h>
#include "ARPACK.h"
#include "MatOp.h"

using Rcpp::as;
using Rcpp::wrap;

// Main function to calculate real, nonsymmetric SVD
// If A = UDV', then A'A = V(D^2)V', AA' = U(D^2)U'
// So if m >= n, we calculate V first, by solving A'A * v = lambda * v
// If m < n, we calculate U first, by solving AA' * u = lambda * u
SEXP do_svds_gen(MatOp *op, SEXP m_scalar_r, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r)
{
BEGIN_RCPP

    // Retrieve parameters
    Rcpp::List params_rcpp(params_list_r);

    // Dimension of A (m by n)
    int m = INTEGER(m_scalar_r)[0];
    int n = INTEGER(n_scalar_r)[0];
    // If m > n, OP = A'A, otherwise OP = AA'
    // The actual dimension of problem is min(m, n)
    int dim = m > n ? n : m;
    // Number of returned singular vectors
    int nu = INTEGER(nu_scalar_r)[0];
    int nv = INTEGER(nv_scalar_r)[0];

    // Begin ARPACK
    //
    // Initial value of ido
    int ido = 0;
    // 'I' means standard eigen value problem, A * x = lambda * x
    char bmat = 'I';
    // Specify selection criteria
    // For SVD we only use "LM"
    char which[3];
    which[0] = 'L';
    which[1] = 'M';
    which[2] = '\0';
    // Number of singular values requested
    int nev = INTEGER(k_scalar_r)[0];
    // Precision
    double tol = as<double>(params_rcpp["tol"]);
    // Residual vector
    double *resid = new double[dim]();
    int ncoef = m > n ? m : n;
    double *tmp = new double[ncoef]();
    for(int i = 0; i < ncoef; i++)
        tmp[i] = sin(i + 0.5);
    if (m > n)
        op->tprod(tmp, resid);
    else
        op->prod(tmp, resid);
    // Related to the algorithm, large ncv results in
    // faster convergence, but with greater memory use
    int ncv = as<int>(params_rcpp["ncv"]);
    
    // Variables to be returned to R
    //
    // Vector of singular values
    Rcpp::NumericVector d_ret(nev);
    // Result list
    Rcpp::List ret;
    
    // Store eigenvectors of OP
    Rcpp::NumericMatrix ev(dim, ncv);
    double *V = ev.begin();
    // Leading dimension of V, required by FORTRAN
    int ldv = dim;
    // Control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1; // ishfts
    iparam[3 - 1] = as<int>(params_rcpp["maxitr"]); // maxitr
    iparam[7 - 1] = 1; // mode
    // Some pointers
    int *ipntr = new int[11]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * dim]();
    int lworkl = ncv * (ncv + 8);
    double *workl = new double[lworkl]();
    // Error flag. 0 means random initialization,
    // otherwise using resid as initial value
    int info = 1;

    saupd(ido, bmat, dim, which,
          nev, tol, resid,
          ncv, V, ldv,
          iparam, ipntr, workd,
          workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        if (m > n)
        {
            // OP = A'A
            // First do tmp = A * x_in
            op->prod(&workd[ipntr[0] - 1], tmp);
            // Then do y_out = A' * tmp
            op->tprod(tmp, &workd[ipntr[1] - 1]);
        } else {
            // OP = AA'
            // First do tmp = A' * x_in
            op->tprod(&workd[ipntr[0] - 1], tmp);
            // Then do y_out = A * tmp
            op->prod(tmp, &workd[ipntr[1] - 1]);
        }
        saupd(ido, bmat, dim, which,
              nev, tol, resid,
              ncv, V, ldv,
              iparam, ipntr, workd,
              workl, lworkl, info);
    }
    delete [] tmp;  
    // info > 0 means warning, < 0 means error
    if (info > 0) dsaupd_warn(info);
    if (info < 0)
    {
        delete [] workl;
        delete [] workd;
        delete [] ipntr;
        delete [] iparam;
        delete [] resid;
        dsaupd_error(info);
    }
    
    // Retrieve results
    //
    // Whether to calculate singular vectors or not.
    bool rvec = (nu > 0) | (nv > 0);
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Vector of singular values
    double *d = d_ret.begin();
    // Used to store results, will use V instead.
    double *Z = V;
    // Leading dimension of Z, required by FORTRAN
    int ldz = dim;
    // Shift
    double sigma = 0;
    // Error information
    int ierr = 0;
    
    // Number of converged eigenvalues
    int nconv = 0;
    // Number of iterations
    int niter = 0;

    // Use seupd() to retrieve results
    seupd(rvec, howmny, d,
          Z, ldz, sigma, bmat,
          dim, which, nev, tol,
          resid, ncv, V, ldv,
          iparam, ipntr, workd, workl,
          lworkl, ierr);
    // Obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    // 'niter' number of iterations
    niter = iparam[9 - 1];

    // Free memory of temp arrays
    delete [] workl;
    delete [] workd;
    delete [] ipntr;
    delete [] iparam;
    delete [] resid;

    // ierr < 0 means error
    if (ierr < 0)  dseupd_error(ierr);

    if (nconv <= 0) ::Rf_error("no converged singular values found");

    // Give a warning if less singular values than requested converged
    if (nconv < nev)
        ::Rf_warning("only %d singular values converged, less than k", nconv);

    // Singular value is the square root of the eigenvalue of OP
    Rcpp::NumericVector s_ret = Rcpp::sqrt(d_ret);
    // ARPACK gives eigenvalues in increasing order.
    // We need decreasing one.
    std::reverse(s_ret.begin(), s_ret.end());
    // If we don't need singular vectors
    if (!rvec)
    {
         ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                  Rcpp::Named("u") = R_NilValue,
                                  Rcpp::Named("v") = R_NilValue,
                                  Rcpp::Named("nconv") = wrap(nconv),
                                  Rcpp::Named("niter") = wrap(niter));
    } else if (m > n) {
        // Change the order of columns of ev, since singular values are reversed
        for(int i = 0; i < nev / 2; i++)
        {
            std::swap_ranges(&ev(0, i), &ev(0, i + 1), &ev(0, nev - i - 1));
        }
        // A = UDV', OP = A'A
        // AV = UD, u = Av / sigma
        // ev contains V, we are now computing u_ret which contains U
        if (nu == 0)
        {
            // Here nv will not be 0
            ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                     Rcpp::Named("u") = R_NilValue,
                                     Rcpp::Named("v") = ev(Rcpp::_, Rcpp::Range(0, nv - 1)),
                                     Rcpp::Named("nconv") = wrap(nconv),
                                     Rcpp::Named("niter") = wrap(niter));
        } else {
            Rcpp::NumericMatrix u_ret(m, nu);
            double *u = u_ret.begin();
            double *v = ev.begin();
            for (int i = 0; i < nu; i++)
            {
                op->prod(v, u);
                for (int j = 0; j < m; j++)
                    u[j] /= s_ret[i];
                u += m;
                v += n;
            }
            if (nv == 0)
                ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                         Rcpp::Named("u") = u_ret,
                                         Rcpp::Named("v") = R_NilValue,
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(niter));
            else
                ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                         Rcpp::Named("u") = u_ret,
                                         Rcpp::Named("v") = ev(Rcpp::_, Rcpp::Range(0, nv - 1)),
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(niter));
        }
    } else {
        // Change the order of columns of ev, since singular values are reversed
        for(int i = 0; i < nev / 2; i++)
        {
            std::swap_ranges(&ev(0, i), &ev(0, i + 1), &ev(0, nev - i - 1));
        }
        // A = UDV', OP = AA'
        // A'U = VD, v = A'u / sigma
        // ev contains U, we are now computing v_ret which contains V
        if (nv == 0)
        {
            // Here nu will not be 0
            ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                     Rcpp::Named("u") = ev(Rcpp::_, Rcpp::Range(0, nu - 1)),
                                     Rcpp::Named("v") = R_NilValue,
                                     Rcpp::Named("nconv") = wrap(nconv),
                                     Rcpp::Named("niter") = wrap(niter));
        } else {
            Rcpp::NumericMatrix v_ret(n, nv);
            double *u = ev.begin();
            double *v = v_ret.begin();
            for (int i = 0; i < nv; i++)
            {
                op->tprod(u, v);
                for (int j = 0; j < n; j++)
                    v[j] /= s_ret[i];
                u += m;
                v += n;
            }
            if (nu == 0)
                ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                         Rcpp::Named("u") = R_NilValue,
                                         Rcpp::Named("v") = v_ret,
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(niter));
            else
                ret = Rcpp::List::create(Rcpp::Named("d") = s_ret,
                                         Rcpp::Named("u") = ev(Rcpp::_, Rcpp::Range(0, nu - 1)),
                                         Rcpp::Named("v") = v_ret,
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(niter));
        }
    }

    return ret;

END_RCPP
}

