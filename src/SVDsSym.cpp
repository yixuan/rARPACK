#include "SVDsSym.h"

SVDsSym::SVDsSym(int n_, int k_, int nu_, int nv_, int ncv_, MatOp *op_,
                 double tol_, int maxitr_) :
    EigsSym(n_, k_, ncv_, op_, "LM", 1L, 'I', tol_, maxitr_),
    nu(nu_), nv(nv_)
{}

Rcpp::List SVDsSym::extract()
{
    // Obtain 'nconv' converged singular values
    int nconv = iparam[5 - 1];
    // 'niter' number of iterations
    int niter = iparam[9 - 1];

    if(nconv <= 0)
    {
        ::Rf_warning("no converged singular values found");
        return returnResult(R_NilValue, R_NilValue, R_NilValue,
                            Rcpp::wrap(nconv), Rcpp::wrap(niter));
    }
    
    if(nconv < nev)
        ::Rf_warning("only %d singular values converged, less than k = %d",
                     nconv, nev);
    
    // Calculate singular values
    eigd.erase(nconv, eigd.length());
    // Sort singular values in decreasing order
    Rcpp::IntegerVector order = sort_with_order<DESCEND>(eigd);
    
    // Copy singular vectors
    Rcpp::RObject U, V;
    Rcpp::NumericMatrix matU, matV;
    nu = nu > nconv ? nconv : nu;
    nv = nv > nconv ? nconv : nv;
    
    if(nu == 0)
    {
        U = R_NilValue;
    } else {
        matU = Rcpp::NumericMatrix(n, nu);
        for(int i = 0; i < nu; i++)
        {
            copy_column(eigV, order[i], matU, i);
        }
        U = matU;
    }
    
    if(nv == 0)
    {
        V = R_NilValue;
    } else {
        matV = Rcpp::NumericMatrix(n, nv);
        for(int i = 0; i < nv; i++)
        {
            copy_column(eigV, order[i], matV, i);
        }
        V = matV;
    }

    return returnResult(eigd, U, V, Rcpp::wrap(nconv), Rcpp::wrap(niter));
}
