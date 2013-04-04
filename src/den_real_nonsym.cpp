#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include "arpack_cpp.h"

// calculate matrix(n, n) * vector(n) product using BLAS
void mat_v_prod(double *mat, double *x, int len, double *y)
{
    char trans = 'N';
    int n = len;
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dgemv)(&trans, &n, &n,
            &alpha, mat, &n,
            x, &one, &zero,
            y, &one);
}


RcppExport SEXP den_real_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP sigmar_scalar_r, SEXP sigmai_scalar_r)
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
    char *which = new char[3];
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
    // vector of real part of eigenvalues
    Rcpp::NumericVector dreal_ret(nev + 1);
    // vector of imag part of eigenvalues
    Rcpp::NumericVector dimag_ret(nev + 1);
    // matrix of real part of eigenvectors
    Rcpp::NumericMatrix vreal_ret(n, ncv);
    // matrix of imag part of eigenvectors
    // determine the dimension of vimag_ret after
    // we get 'nconv' converged eigenvalues
    Rcpp::NumericMatrix vimag_ret;
    // result list
    Rcpp::List ret;
    
    // store final results of eigenvectors
    // double *V = new double[n * ncv]();
    double *V = vreal_ret.begin();
    // leading dimension of V, required by FORTRAN
    int ldv = n;
    // control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1; // ishfts
    iparam[3 - 1] = INTEGER(maxitr_scalar_r)[0]; // maxitr
    iparam[7 - 1] = 1; // mode
    // some pointers
    int *ipntr = new int[14]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * n]();
    int lworkl = 3 * ncv * ncv + 6 * ncv;
    double *workl = new double[lworkl]();
    // error flag, initialized to 0
    int info = 0;

    naupp(ido, bmat, n, which, nev,
            tol, resid, ncv, V,
            ldv, iparam, ipntr, workd,
            workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        mat_v_prod(REAL(A_mat_r), &workd[ipntr[0] - 1], n,
                &workd[ipntr[1] - 1]);
        naupp(ido, bmat, n, which, nev,
            tol, resid, ncv, V,
            ldv, iparam, ipntr, workd,
            workl, lworkl, info);
    }
    
    // retrieve results
    //
    // whether to calculate eigenvectors or not.
    bool rvec = true;
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char HowMny = 'A';
    // real part of eigenvalues
    double *dr = dreal_ret.begin();
    // imaginary part of eigenvalues
    double *di = dimag_ret.begin();
    // used to store results, will use V instead.
    double *Z = V;
    // leading dimension of Z, required by FORTRAN
    int ldz = n;
    // shift
    double sigmar = REAL(sigmar_scalar_r)[0];
    double sigmai = REAL(sigmai_scalar_r)[0];
    // working space
    double *workv = new double[3 * ncv]();
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
        neupp(rvec, HowMny, dr,
                di, Z, ldz, sigmar,
                sigmai, workv, bmat, n,
                which, nev, tol, resid,
                ncv, V, ldv, iparam,
                ipntr, workd, workl,
                lworkl, ierr);
        if (ierr < 0)
        {
            ::Rf_error("Error in dneupd subroutine of ARPACK,"
                       "with code %d", ierr);
        } else {
            // obtain 'nconv' converged eigenvalues
            nconv = iparam[5 - 1];
            Rcpp::Range range = Rcpp::Range(0, nconv - 1);

            // if all eigenvalues are real
            // equivalent R code: if (all(abs(dimag[1:nconv] < 1e-30)))
            if (Rcpp::is_true(Rcpp::all(Rcpp::abs(dimag_ret) < 1e-30)))
            {
                // so we don't need to allocate 'vimag_ret'
                // for imaginary part of the eigenvectors
                dreal_ret.erase(nconv, dreal_ret.length() - 1);
                ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                         Rcpp::Named("values") = dreal_ret,
                                         Rcpp::Named("vectors") = vreal_ret(Rcpp::_, range));
            } else {
                vimag_ret = Rcpp::NumericMatrix(n, nconv);
                Rcpp::ComplexVector cmpvalues_ret(nconv);
                Rcpp::ComplexMatrix cmpvectors_ret(n, nconv);
                // obtain the imaginary part of the eigenvectors
                bool first = true;
                for (int i = 0; i < nconv; i++)
                {
                    cmpvalues_ret[i].r = dreal_ret[i];
                    cmpvalues_ret[i].i = dimag_ret[i];
                    if (first && fabs(dimag_ret[i]) > 1e-30)
                    {
                        vimag_ret(Rcpp::_, i) = vreal_ret(Rcpp::_, i + 1);
                        vimag_ret(Rcpp::_, i + 1) = -vreal_ret(Rcpp::_, i + 1);
                        vreal_ret(Rcpp::_, i + 1) = vreal_ret(Rcpp::_, i);
                        first = false;
                    } else {
                        first = true;
                    }
                }
                // assign values to complex eigenvectors
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < nconv; j++)
                    {
                        cmpvectors_ret(i, j).r = vreal_ret(i, j);
                        cmpvectors_ret(i, j).i = vimag_ret(i, j);
                    }
                }
                ret = Rcpp::List::create(Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                         Rcpp::Named("values") = cmpvalues_ret,
                                         Rcpp::Named("vectors") = cmpvectors_ret);
            }
        }
    }

    delete [] which;
    delete [] workv;
    delete [] workl;
    delete [] workd;
    delete [] ipntr;
    delete [] iparam;
    delete [] resid;

    return ret;

END_RCPP
}

