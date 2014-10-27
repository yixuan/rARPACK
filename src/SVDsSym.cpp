#include "SVDsSym.h"

SVDsSym::SVDsSym(int n_, int k_, int nu_, int nv_, int ncv_, MatOp *op_,
                 double tol_, int maxitr_) :
    EigsSym(n_, k_, ncv_, op_, "LM", 1L, 'I', tol_, maxitr_),
    nu(nu_), nv(nv_)
{}

void SVDsSym::compute()
{
    bool rvec = (nu > 0) || (nv > 0);
    EigsSym::compute(rvec);
}

Rcpp::RObject SVDsSym::extractEigenvectors(int num)
{
    if(num <= 0)
        return R_NilValue;
    
    int nconv = iparam[5 - 1];
    if(num > nconv)  num = nconv;
    
    Rcpp::NumericMatrix mat(n, num);
    for(int i = 0; i < num; i++)
    {
        copy_column(eigV, dind[i], mat, i);
    }
    
    return mat;
}

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
    dind = sort_with_order<ABSDESCEND>(eigd);
    dval = eigd;
    
    // Copy singular vectors
    Rcpp::RObject U = extractEigenvectors(nu);
    Rcpp::RObject V = extractEigenvectors(nv);
    
    // We need to make sure that singular values are nonnegative,
    // so move the sign to matV.
    for(int i = 0; i < nconv; i++)
    {
        if(dval[i] < 0)
        {
            dval[i] = -dval[i];
            if(i < nv)
            {
                SEXP Vdata = V;
                Rcpp::NumericMatrix matV(Vdata);
                matV(Rcpp::_, i) = -matV(Rcpp::_, i);
            }
        }
    }

    return returnResult(dval, U, V, Rcpp::wrap(nconv), Rcpp::wrap(niter));
}
