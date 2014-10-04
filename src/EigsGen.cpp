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
    int truenconv = nconv > nev ? nev : nconv;

    // equivalent R code: if(all(abs(dimag[1:nconv] < 1e-17)))
    if(Rcpp::is_true(Rcpp::all(Rcpp::abs(eigdi) < 1e-17)))
    {
        // v.erase(start, end) removes v[start <= i < end]
        eigdr.erase(truenconv, eigdr.length());
        // Make sure eigenvalues are in decreasing order.
        Rcpp::IntegerVector order = sort_with_order(eigdr);
        
        if(!retvec)
        {
            ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
            return ret;
        }
        
        // The matrix to be returned
        Rcpp::NumericMatrix retV(n, truenconv);
        for(int i = 0; i < truenconv; i++)
        {
            copy_column(eigV, order[i], retV, i);
        }

        ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                 Rcpp::Named("vectors") = retV,
                                 Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
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
        // To make sure the eigenvalues are in proper order
        sort_with_order(cmpeigd);
        if(nconv > truenconv)  cmpeigd.erase(truenconv, nconv);
        ret = Rcpp::List::create(Rcpp::Named("values") = cmpeigd,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
        return ret;
    }
    
    Rcpp::ComplexMatrix retV(n, truenconv);
    Rcpp::ComplexVector retd = Rcpp::clone(cmpeigd);
    bool *img_zero = new bool[nconv + 1];  // add 1 to be safe
    for(int i = 0; i < nconv; i++)
    {
        if(fabs(eigdi[i]) < 1e-16)
        {
            img_zero[i] = true;
        } else {
            img_zero[i] = false;
            i++;
            img_zero[i] = false;
            retd[i].r = 0;
            retd[i].i = 0;
        }
    }
    Rcpp::IntegerVector order = sort_with_order(retd);
    int *order_pntr = order.begin();
    for(int i = 0; i < truenconv; i++, order_pntr++)
    {
        int ind = *order_pntr;
        if(img_zero[ind])
        {
            retd[i] = cmpeigd[ind];
            copy_column(eigV, ind, eigV, ind, 0, retV, i);
        } else {
            retd[i] = cmpeigd[ind];
            copy_column(eigV, ind, eigV, ind + 1, 1, retV, i);
            i++;
            if(i < truenconv)
            {
                retd[i].r = cmpeigd[ind].r;
                retd[i].i = -cmpeigd[ind].i;
                copy_column(eigV, ind, eigV, ind + 1, -1, retV, i);
            }
        }
    }
    delete [] img_zero;
    if(nconv > truenconv)  retd.erase(truenconv, nconv);
    
    ret = Rcpp::List::create(Rcpp::Named("values") = retd,
                             Rcpp::Named("vectors") = retV,
                             Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                             Rcpp::Named("niter") = Rcpp::wrap(niter));

    return ret;
}
