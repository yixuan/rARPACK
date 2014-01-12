#include "do_eigs.h"

using Rcpp::as;
using Rcpp::wrap;

SEXP do_eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                 SEXP params_list_r,
                 Mvfun mat_v_prod, void *data)
{
BEGIN_RCPP

    // Retrieve parameters
    Rcpp::List params_rcpp(params_list_r);

    // Begin ARPACK
    //
    // Initial value of ido
    int ido = 0;
    // 'I' means standard eigen value problem, A * x = lambda * x
    char bmat = 'I';
    // dimension of A (n by n)
    int n = INTEGER(n_scalar_r)[0];
    // Specify selection criteria
    // "LM": largest magnitude
    // "SM": smallest magnitude
    // "LR", "LI": largest real/imaginary part
    // "SR", "SI": smallest real/imaginary part
    Rcpp::CharacterVector which_rcpp = params_rcpp["which"];
    char which[3];
    which[0] = which_rcpp[0][0];
    which[1] = which_rcpp[0][1];
    which[2] = '\0';
    // Number of eigenvalues requested
    int nev = INTEGER(k_scalar_r)[0];
    // Precision
    double tol = as<double>(params_rcpp["tol"]);
    // Residual vector
    double *resid = new double[n]();
    double *initcoef = new double[n]();
    initcoef[0] = initcoef[n - 1] = 0.5;
    mat_v_prod(A_mat_r, initcoef, resid, n, n, data);
    delete [] initcoef;
    // Related to the algorithm, large ncv results in
    // faster convergence, but with greater memory use
    int ncv = as<int>(params_rcpp["ncv"]);
    
    // Variables to be returned to R
    //
    // Vector of eigenvalues
    Rcpp::NumericVector d_ret(nev);
    // Matrix of eigenvectors
    Rcpp::NumericMatrix v_ret(n, ncv);
    // Result list
    Rcpp::List ret;
    
    // Store final results of eigenvectors
    // double *V = new double[n * ncv]();
    double *V = v_ret.begin();
    // Leading dimension of V, required by FORTRAN
    int ldv = n;
    // Control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1; // ishfts
    iparam[3 - 1] = as<int>(params_rcpp["maxitr"]); // maxitr
    iparam[7 - 1] = as<int>(params_rcpp["workmode"]); // mode
    // Some pointers
    int *ipntr = new int[11]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * n]();
    int lworkl = ncv * (ncv + 8);
    double *workl = new double[lworkl]();
    // Error flag. 0 means random initialization,
    // otherwise using resid as initial value
    int info = 1;

    saupd(ido, bmat, n, which,
          nev, tol, resid,
          ncv, V, ldv,
          iparam, ipntr, workd,
          workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        mat_v_prod(A_mat_r, &workd[ipntr[0] - 1],
                   &workd[ipntr[1] - 1], n, n, data);
        saupd(ido, bmat, n, which,
              nev, tol, resid,
              ncv, V, ldv,
              iparam, ipntr, workd,
              workl, lworkl, info);
    }
    
    // info > 0 means warning, < 0 means error
    if(info > 0) dsaupd_warn(info);
    if(info < 0)
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
    // Whether to calculate eigenvectors or not.
    bool rvec = as<bool>(params_rcpp["retvec"]);
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Vector of eigenvalues
    double *d = d_ret.begin();
    // Used to store results, will use V instead.
    double *Z = V;
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    // Shift
    double sigma = as<double>(params_rcpp["sigma"]);
    // Error information
    int ierr = 0;
    
    // Number of converged eigenvalues
    int nconv = 0;

    // Use seupd() to retrieve results
    seupd(rvec, howmny, d,
          Z, ldz, sigma, bmat,
          n, which, nev, tol,
          resid, ncv, V, ldv,
          iparam, ipntr, workd, workl,
          lworkl, ierr);

    // ierr < 0 means error
    if (ierr < 0)
    {
        delete [] workl;
        delete [] workd;
        delete [] ipntr;
        delete [] iparam;
        delete [] resid;
        dseupd_error(ierr);
    }
    
    // Obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    if(nconv <= 0)
    {
         ::Rf_warning("no converged eigenvalues found");
         ret = Rcpp::List::create(Rcpp::Named("values") = R_NilValue,
                                  Rcpp::Named("vectors") = R_NilValue,
                                  Rcpp::Named("nconv") = wrap(nconv),
                                  Rcpp::Named("niter") = wrap(iparam[9 - 1]));
    } else {
        if (nconv < nev)
            ::Rf_warning("only %d eigenvalues converged, less than k", nconv);
        // v.erase(start, end) removes v[start <= i < end]
        d_ret.erase(nconv, d_ret.length());
        // ARPACK gives eigenvalues in increasing order.
        // We need decreasing one.
        std::reverse(d_ret.begin(), d_ret.end());
        if(rvec)
        {
            // Also change the order of columns of v_ret
            for(int i = 0; i < nconv / 2; i++)
            {
                std::swap_ranges(&v_ret(0, i), &v_ret(0, i + 1),
                                 &v_ret(0, nconv - i - 1));
            }
            Rcpp::Range range = Rcpp::Range(0, nconv - 1);
            ret = Rcpp::List::create(Rcpp::Named("values") = d_ret,
                                     Rcpp::Named("vectors") = v_ret(Rcpp::_, range),
                                     Rcpp::Named("nconv") = wrap(nconv),
                                     Rcpp::Named("niter") = wrap(iparam[9 - 1]));
        } else {
            ret = Rcpp::List::create(Rcpp::Named("values") = d_ret,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("nconv") = wrap(nconv),
                                     Rcpp::Named("niter") = wrap(iparam[9 - 1]));
        }
    }

    delete [] workl;
    delete [] workd;
    delete [] ipntr;
    delete [] iparam;
    delete [] resid;

    return ret;

END_RCPP
}

