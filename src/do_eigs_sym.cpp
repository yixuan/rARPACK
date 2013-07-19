#include "do_eigs.h"

SEXP do_eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r, SEXP sigma_scalar_r,
        Mvfun mat_v_prod, void *data)
{
BEGIN_RCPP

    // begin ARPACK
    //
    // initial value of ido
    int ido = 0;
    // 'I' means standard eigen value problem, A * x = lambda * x
    char bmat = 'I';
    // dimension of A (n by n)
    int n = INTEGER(n_scalar_r)[0];
    // specify selection criteria
    // "LM": largest magnitude
    // "SM": smallest magnitude
    // "LR", "LI": largest real/imaginary part
    // "SR", "SI": smallest real/imaginary part
    char which[3];
    which[0] = CHAR(STRING_ELT(which_string_r, 0))[0];
    which[1] = CHAR(STRING_ELT(which_string_r, 0))[1];
    which[2] = '\0';
    // number of eigenvalues requested
    int nev = INTEGER(k_scalar_r)[0];
    // precision
    double tol = REAL(tol_scalar_r)[0];
    // residual vector
    double *resid = new double[n]();
    // related to the algorithm, large ncv results in
    // faster convergence, but with greater memory use
    int ncv = INTEGER(ncv_scalar_r)[0];
    
    // variables to be returned to R
    //
    // vector of eigenvalues
    Rcpp::NumericVector d_ret(nev + 1);
    // matrix of eigenvectors
    Rcpp::NumericMatrix v_ret(n, ncv);
    // result list
    Rcpp::List ret;
    
    // store final results of eigenvectors
    // double *V = new double[n * ncv]();
    double *V = v_ret.begin();
    // leading dimension of V, required by FORTRAN
    int ldv = n;
    // control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1; // ishfts
    iparam[3 - 1] = INTEGER(maxitr_scalar_r)[0]; // maxitr
    iparam[7 - 1] = 1; // mode
    // some pointers
    int *ipntr = new int[11]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * n]();
    int lworkl = ncv * (ncv + 8);
    double *workl = new double[lworkl]();
    // error flag, initialized to 0
    int info = 0;

    saupp(ido, bmat, n, which, nev,
            tol, resid, ncv, V,
            ldv, iparam, ipntr, workd,
            workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        mat_v_prod(A_mat_r, &workd[ipntr[0] - 1],
                   &workd[ipntr[1] - 1], n, data);
        saupp(ido, bmat, n, which, nev,
            tol, resid, ncv, V,
            ldv, iparam, ipntr, workd,
            workl, lworkl, info);
    }
    
    // retrieve results
    //
    // whether to calculate eigenvectors or not.
    bool rvec = (bool) LOGICAL(retvec_logical_r)[0];
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char HowMny = 'A';
    // vector of eigenvalues
    double *d = d_ret.begin();
    // used to store results, will use V instead.
    double *Z = V;
    // leading dimension of Z, required by FORTRAN
    int ldz = n;
    // shift
    double sigma = REAL(sigma_scalar_r)[0];
    // error information
    int ierr = 0;
    
    // number of converged eigenvalues
    int nconv = 0;

    // info < 0 means error occurs
    if (info < 0)
    {
        ::Rf_error("Error in dnaupd subroutine of ARPACK, with code %d",
                   info);
    } else {
        // use neupp() to retrieve results
        seupp(rvec, HowMny, d, Z, ldz, sigma,
                bmat, n,
                which, nev, tol, resid,
                ncv, V, ldv, iparam,
                ipntr, workd, workl,
                lworkl, ierr);
        if (ierr < 0)
        {
            ::Rf_error("Error in dseupd subroutine of ARPACK,"
                       "with code %d", ierr);
        } else {
            // obtain 'nconv' converged eigenvalues
            nconv = iparam[5 - 1];
            Rcpp::Range range = Rcpp::Range(0, nconv - 1);

            d_ret.erase(nconv, d_ret.length() - 1);
            if(rvec)
            {
                ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("values") = d_ret,
                                     Rcpp::Named("vectors") = v_ret(Rcpp::_, range));
            } else {
                ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("values") = d_ret,
                                     Rcpp::Named("vectors") = R_NilValue);
            }
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

