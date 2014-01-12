#include "do_eigs.h"

using Rcpp::as;
using Rcpp::wrap;

SEXP do_eigs_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
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
    // Dimension of A (n by n)
    int n = INTEGER(n_scalar_r)[0];
    // specify selection criteria
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
    // Vector of real part of eigenvalues
    Rcpp::NumericVector dreal_ret(nev + 1);
    // Vector of imag part of eigenvalues
    Rcpp::NumericVector dimag_ret(nev + 1);
    // Matrix of real part of eigenvectors
    Rcpp::NumericMatrix vreal_ret(n, ncv);
    // Result list
    Rcpp::List ret;
    
    // Store final results of eigenvectors
    // double *V = new double[n * ncv]();
    double *V = vreal_ret.begin();
    // Leading dimension of V, required by FORTRAN
    int ldv = n;
    // Control parameters
    int *iparam = new int[11]();
    iparam[1 - 1] = 1; // ishfts
    iparam[3 - 1] = as<int>(params_rcpp["maxitr"]); // maxitr
    iparam[7 - 1] = as<int>(params_rcpp["workmode"]); // mode
    // Some pointers
    int *ipntr = new int[14]();
    /* workd has 3 columns.
     * ipntr[2] - 1 ==> first column to store B * X,
     * ipntr[1] - 1 ==> second to store Y,
     * ipntr[0] - 1 ==> third to store X. */
    double *workd = new double[3 * n]();
    int lworkl = 3 * ncv * ncv + 6 * ncv;
    double *workl = new double[lworkl]();
    // Error flag. 0 means random initialization,
    // otherwise using resid as initial value
    int info = 1;

    naupd(ido, bmat, n, which,
          nev, tol, resid,
          ncv, V, ldv,
          iparam, ipntr, workd,
          workl, lworkl, info);
    // ido == -1 or ido == 1 means more iterations needed
    while (ido == -1 || ido == 1)
    {
        mat_v_prod(A_mat_r, &workd[ipntr[0] - 1],
                   &workd[ipntr[1] - 1], n, n, data);
        naupd(ido, bmat, n, which,
              nev, tol, resid,
              ncv, V, ldv,
              iparam, ipntr, workd,
              workl, lworkl, info);
    }
    
    // info > 0 means warning, < 0 means error
    if(info > 0) dnaupd_warn(info);
    if(info < 0)
    {
        delete [] workl;
        delete [] workd;
        delete [] ipntr;
        delete [] iparam;
        delete [] resid;
        dnaupd_error(info);
    }
    
    // Retrieve results
    //
    // Whether to calculate eigenvectors or not.
    bool rvec = as<bool>(params_rcpp["retvec"]);
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Real part of eigenvalues
    double *dr = dreal_ret.begin();
    // Imaginary part of eigenvalues
    double *di = dimag_ret.begin();
    // Used to store results, will use V instead.
    double *Z = V;
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    // Shift
    double sigmar = as<double>(params_rcpp["sigmar"]);
    double sigmai = as<double>(params_rcpp["sigmai"]);
    // Working space
    double *workv = new double[3 * ncv]();
    // Error information
    int ierr = 0;
    
    // Number of converged eigenvalues
    int nconv = 0;

    // Use neupd() to retrieve results
    neupd(rvec, howmny, dr, di,
          Z, ldz, sigmar, sigmai, workv,
          bmat, n, which, nev, tol,
          resid, ncv, V, ldv, iparam,
          ipntr, workd, workl, lworkl, ierr);
    // ierr > 0 means warning, < 0 means error
    if (ierr > 0) dneupd_warn(ierr);
    if (ierr < 0)
    {
        delete [] workv;
        delete [] workl;
        delete [] workd;
        delete [] ipntr;
        delete [] iparam;
        delete [] resid;
        dneupd_error(ierr);
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
        // Sometimes there are nconv = nev + 1 converged eigenvalues,
        // mainly due to pairs of complex eigenvalues.
        // We will truncate at nev
        if (nconv > nev)  nconv = nev;
        // If all eigenvalues are real
        // equivalent R code: if (all(abs(dimag[1:nconv] < 1e-17)))
        if (Rcpp::is_true(Rcpp::all(Rcpp::abs(dimag_ret) < 1e-17)))
        {
            // v.erase(start, end) removes v[start <= i < end]
            dreal_ret.erase(nconv, dreal_ret.length());
            if(rvec)
            {
                Rcpp::Range range = Rcpp::Range(0, nconv - 1);
                ret = Rcpp::List::create(Rcpp::Named("values") = dreal_ret,
                                         Rcpp::Named("vectors") = vreal_ret(Rcpp::_, range),
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(iparam[9 - 1]));
            } else {
                ret = Rcpp::List::create(Rcpp::Named("values") = dreal_ret,
                                         Rcpp::Named("vectors") = R_NilValue,
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(iparam[9 - 1]));
            }        
        } else {
            Rcpp::ComplexVector cmpvalues_ret(nconv);
            int i;
            for (i = 0; i < nconv; i++)
            {
                cmpvalues_ret[i].r = dreal_ret[i];
                cmpvalues_ret[i].i = dimag_ret[i];
            }
            if(!rvec)
            {
                // It seems that when rvec is false, ARPACK will return
                // the eigenvalues in ascending order. So here we reverse
                // the vector
                std::reverse(cmpvalues_ret.begin(), cmpvalues_ret.end());
                ret = Rcpp::List::create(Rcpp::Named("values") = cmpvalues_ret,
                                         Rcpp::Named("vectors") = R_NilValue,
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(iparam[9 - 1]));
            } else {
                Rcpp::ComplexMatrix cmpvectors_ret(n, nconv);
                // Obtain the real and imaginary part of the eigenvectors
                bool first = true;
                int j;
                for (i = 0; i < nconv; i++)
                {
                    if (fabs(dimag_ret[i]) > 1e-30)
                    {
                        if (first)
                        {
                            for (j = 0; j < n; j++)
                            {
                                cmpvectors_ret(j, i).r = vreal_ret(j, i);
                                cmpvectors_ret(j, i).i = vreal_ret(j, i + 1);
                                cmpvectors_ret(j, i + 1).r = vreal_ret(j, i);
                                cmpvectors_ret(j, i + 1).i = -vreal_ret(j, i + 1);
                            }
                            first = false;
                        } else {
                            first = true;
                        }
                    } else {
                        for (j = 0; j < n; j++)
                        {
                            cmpvectors_ret(j, i).r = vreal_ret(j, i);
                            cmpvectors_ret(j, i).i = 0;
                        }
                        first = true;
                    }
                }
                ret = Rcpp::List::create(Rcpp::Named("values") = cmpvalues_ret,
                                         Rcpp::Named("vectors") = cmpvectors_ret,
                                         Rcpp::Named("nconv") = wrap(nconv),
                                         Rcpp::Named("niter") = wrap(iparam[9 - 1]));
            }
            
        }
    }
    

    delete [] workv;
    delete [] workl;
    delete [] workd;
    delete [] ipntr;
    delete [] iparam;
    delete [] resid;

    return ret;

END_RCPP
}

