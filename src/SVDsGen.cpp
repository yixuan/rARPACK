#include "SVDsGen.h"

using Rcpp::as;
using Rcpp::wrap;

Rcpp::RObject SVDsGenTall::computeU(int num)
{
    if(num <= 0)
        return R_NilValue;
    
    int nconv = iparam[5 - 1];
    if(num > nconv)  num = nconv;
    
    // AV = UD, u = Av / d
    int nrowU = nrow > ncol ? nrow : ncol;
    Rcpp::NumericMatrix matU(nrowU, num);
    double *u = matU.begin(), *v;
    for(int i = 0; i < num; i++)
    {
        v = &eigV(0, dind[i]);
        // Calculating A * v
        matProd(v, u);
        // Calculating Av / d
        matU(Rcpp::_, i) = matU(Rcpp::_, i) / dval[i];
        u += nrowU;
    }
    
    return matU;
}

Rcpp::List SVDsGenTall::extract()
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
    dind = sort_with_order<DESCEND>(eigd);
    dval = Rcpp::sqrt(eigd);
    
    // Singular vectors
    // A = UDV', OP = A'A
    // AV = UD, u = Av / d
    Rcpp::RObject V = extractEigenvectors(nv);
    Rcpp::RObject U = computeU(nu);
    
    return returnResult(dval, U, V, Rcpp::wrap(nconv), Rcpp::wrap(niter));
}



Rcpp::List SVDsGenWide::extract()
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
    dind = sort_with_order<DESCEND>(eigd);
    dval = Rcpp::sqrt(eigd);
    
    // Singular vectors
    // A = UDV', OP = AA'
    // A'U = VD, v = A'u / d
    Rcpp::RObject U = extractEigenvectors(nu);
    Rcpp::RObject V = computeV(nv);
    
    return returnResult(dval, U, V, Rcpp::wrap(nconv), Rcpp::wrap(niter));
}
