#include "do_eigs.h"

// Warning and error information
// See ARPACK/dsaupd.f for details
static void dsaupd_warn_nonsym(int info)
{
    switch(info)
    {
        case 1:
            ::Rf_warning("ARPACK/dsaupd: maximum number of iterations taken");
            break;
        case 2:
            break;
        case 3:
            ::Rf_warning("ARPACK/dsaupd: no shifts could be applied, try to increase ncv");
            break;
    }
}

static void dsaupd_error_nonsym(int info)
{
    switch(info)
    {
        case -1:
            ::Rf_error("ARPACK/dsaupd: n must be positive");
            break;
        case -2:
            ::Rf_error("ARPACK/dsaupd: k must be positive");
            break;
        case -3:
            ::Rf_error("ARPACK/dsaupd: k < ncv <= n");
            break;
        case -4:
            ::Rf_error("ARPACK/dsaupd: maxitr must be positive");
            break;
        case -5:
            ::Rf_error("ARPACK/dsaupd: which must be one of 'LM', 'SM', 'LA', 'SA', 'BE'");
            break;
        case -6:
            ::Rf_error("ARPACK/dsaupd: error code %d", info);
            break;
        case -7:
            ::Rf_error("ARPACK/dsaupd: length of private work array WORKL is not sufficient");
            break;
        case -8:
            ::Rf_error("ARPACK/dsaupd: error return from trid. eigenvalue calculation\n"
                       "informational error from LAPACK routine dsteqr");
            break;
        case -9:
            ::Rf_error("ARPACK/dsaupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
        case -13:
            ::Rf_error("ARPACK/dsaupd: error code %d", info);
            break;
        case -9999:
            ::Rf_error("ARPACK/dsaupd: couldn't build an Arnoldi factorization");
            break;
        default:
            ::Rf_error("ARPACK/dsaupd: error code %d", info);
            break;
    }
}

// Error information
// See ARPACK/dseupd.f for details
static void dseupd_error_nonsym(int info)
{
    switch(info)
    {
        case -1:
            ::Rf_error("ARPACK/dseupd: n must be positive");
            break;
        case -2:
            ::Rf_error("ARPACK/dseupd: k must be positive");
            break;
        case -3:
            ::Rf_error("ARPACK/dseupd: k < ncv <= n");
            break;
        case -5:
            ::Rf_error("ARPACK/dseupd: which must be one of 'LM', 'SM', 'LA', 'SA', 'BE'");
            break;
        case -6:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
            break;
        case -7:
            ::Rf_error("ARPACK/dseupd: length of private work WORKL array is not sufficient");
            break;
        case -8:
            ::Rf_error("ARPACK/dseupd: error return from trid. eigenvalue calculation\n"
                       "informational error from LAPACK routine dsteqr");
            break;
        case -9:
            ::Rf_error("ARPACK/dseupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
            break;
        case -14:
            ::Rf_error("ARPACK/dseupd: DSAUPD did not find any eigenvalues to sufficient accuracy");
            break;
        case -15:
        case -16:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
        case -17:
            ::Rf_error("ARPACK/dseupd: DSEUPD got a different count of the number of converged Ritz values than DSAUPD got");
            break;
        default:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
            break;
    }
}

SEXP do_eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r, SEXP sigma_scalar_r,
        SEXP workmode_scalar_r,
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
    Rcpp::NumericVector d_ret(nev);
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
    iparam[7 - 1] = INTEGER(workmode_scalar_r)[0]; // mode
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
    
    // info > 0 means warning, < 0 means error
    if(info > 0) dsaupd_warn_nonsym(info);
    if(info < 0)
    {
        delete [] workl;
        delete [] workd;
        delete [] ipntr;
        delete [] iparam;
        delete [] resid;
        dsaupd_error_nonsym(info);
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

    // use neupp() to retrieve results
    seupp(rvec, HowMny, d, Z, ldz, sigma,
            bmat, n,
            which, nev, tol, resid,
            ncv, V, ldv, iparam,
            ipntr, workd, workl,
            lworkl, ierr);

    // ierr < 0 means error
    if (ierr < 0)
    {
        delete [] workl;
        delete [] workd;
        delete [] ipntr;
        delete [] iparam;
        delete [] resid;
        dseupd_error_nonsym(ierr);
    }
    
    // obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    if(nconv <= 0)
    {
         ::Rf_warning("no converged eigenvalues found");
         ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                  Rcpp::Named("values") = R_NilValue,
                                  Rcpp::Named("vectors") = R_NilValue,
                                  Rcpp::Named("niter") = Rcpp::wrap(iparam[9 - 1]));
    } else {
        // v.erase(start, end) removes v[start <= i < end]
        d_ret.erase(nconv, d_ret.length());
        if(rvec)
        {
            Rcpp::Range range = Rcpp::Range(0, nconv - 1);
            ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("values") = d_ret,
                                     Rcpp::Named("vectors") = v_ret(Rcpp::_, range),
                                     Rcpp::Named("niter") = Rcpp::wrap(iparam[9 - 1]));
        } else {
            ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("values") = d_ret,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("niter") = Rcpp::wrap(iparam[9 - 1]));
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

