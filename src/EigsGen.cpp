#include "EigsGen.h"

using std::string;

EigsGen::EigsGen(int n_, int nev_, int ncv_, MatOp *op_,
                 const string & which_, int workmode_,
                 char bmat_, double tol_, int maxitr_) :
    Eigs(n_, nev_, ncv_, op_, which_, workmode_,
         bmat_, tol_, maxitr_),
    eigV(n, ncv), eigdr(nev + 1), eigdi(nev + 1)
{
    lworkl = 3 * ncv * ncv + 6 * ncv;
    workl = new double[lworkl]();
    workv = new double[3 * ncv]();
}


EigsGen::~EigsGen()
{
    delete [] workv;
    delete [] workl;
}

void EigsGen::error(int stage, int errorcode)
{
    if (stage == 1) // dnaupd
    {
        dnaupd_error(errorcode);
    } else { // dneupd
        dneupd_error(errorcode);
    }
}

void EigsGen::warning(int stage, int errorcode)
{
    if (stage == 1) // dnaupd
    {
        dnaupd_warn(errorcode);
    } else { // dneupd
        dneupd_warn(errorcode);
    }
}

void EigsGen::aupd()
{
    naupd(ido, bmat, n, which.c_str(),
          nev, tol, resid,
          ncv, eigV.begin(), n,
          iparam, ipntr, workd,
          workl, lworkl, info);
}

void EigsGen::eupd()
{
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Used to store results, will use V instead.
    double *Z = eigV.begin();
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    
    // Use neupd() to retrieve results
    neupd(retvec, howmny, eigdr.begin(), eigdi.begin(),
          Z, ldz, op->getsigmar(), op->getsigmai(), workv,
          bmat, n, which.c_str(), nev, tol,
          resid, ncv, eigV.begin(), n,
          iparam, ipntr, workd, workl,
          lworkl, ierr);
}

Rcpp::List EigsGen::extract()
{
    // Result list
    Rcpp::List ret;
    // Obtain 'nconv' converged eigenvalues
    int nconv = iparam[5 - 1];
    // 'niter' number of iterations
    int niter = iparam[9 - 1];

    /*********************************************/
    //
    // Case 1: If there is no converged eigenvalue 
    //
    /*********************************************/
    if(nconv <= 0)
    {
        ::Rf_warning("no converged eigenvalues found");
        ret = Rcpp::List::create(Rcpp::Named("values") = R_NilValue,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
        return ret;
    }

    /*************************************/
    //
    // Case 2: If all eigenvalues are real
    //
    /*************************************/
    if(nconv < nev)
        ::Rf_warning("only %d eigenvalues converged, less than k", nconv);
    // Sometimes there are nconv = nev + 1 converged eigenvalues,
    // mainly due to pairs of complex eigenvalues.
    // We will truncate at nev
    if(nconv > nev)  nconv = nev;

    // equivalent R code: if(all(abs(dimag[1:nconv] < 1e-17)))
    if(Rcpp::is_true(Rcpp::all(Rcpp::abs(eigdi) < 1e-17)))
    {
        // v.erase(start, end) removes v[start <= i < end]
        eigdr.erase(nconv, eigdr.length());
        // Make sure eigenvalues are in decreasing order.
        Rcpp::IntegerVector order = sort_with_order(eigdr);
        
        if(!retvec)
        {
            ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
            return ret;
        }
        
        // The matrix to be returned
        Rcpp::NumericMatrix retV(n, nconv);
        for(int i = 0; i < nconv; i++)
        {
            copy_column(eigV, order[i], retV, i);
        }

        ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                 Rcpp::Named("vectors") = retV,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));

        return ret;
    }

    /*************************************/
    //
    // Case 3: Complex eigenvalues
    //
    /*************************************/
    Rcpp::ComplexVector cmpeigd(nconv);
    for(int i = 0; i < nconv; i++)
    {
        cmpeigd[i].r = eigdr[i];
        cmpeigd[i].i = eigdi[i];
    }
    if(!retvec)
    {
        // It seems that when rvec is false, ARPACK will return
        // the eigenvalues in ascending order. So here we reverse
        // the vector
        std::reverse(cmpeigd.begin(), cmpeigd.end());
        ret = Rcpp::List::create(Rcpp::Named("values") = cmpeigd,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
    } else {
        Rcpp::ComplexMatrix cmpeigV(n, nconv);
        // Obtain the real and imaginary part of the eigenvectors
        bool first = true;
        for (int i = 0; i < nconv; i++)
        {
            if (fabs(eigdi[i]) > 1e-30)
            {
                if (first)
                {
                    for (int j = 0; j < n; j++)
                    {
                        cmpeigV(j, i).r = eigV(j, i);
                        cmpeigV(j, i).i = eigV(j, i + 1);
                        if (i + 1 < nconv)
                            cmpeigV(j, i + 1).r = eigV(j, i);
                        if (i + 1 < nconv)
                            cmpeigV(j, i + 1).i = -eigV(j, i + 1);
                    }
                    first = false;
                } else {
                    first = true;
                }
            } else {
                for (int j = 0; j < n; j++)
                {
                    cmpeigV(j, i).r = eigV(j, i);
                    cmpeigV(j, i).i = 0;
                }
                first = true;
            }
        }
        ret = Rcpp::List::create(Rcpp::Named("values") = cmpeigd,
                                 Rcpp::Named("vectors") = cmpeigV,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
    }

    return ret;
}
