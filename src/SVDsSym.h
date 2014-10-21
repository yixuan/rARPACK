#ifndef SVDSSYM_H
#define SVDSSYM_H

#include "EigsSym.h"

class SVDsSym: public EigsSym
{
private:
    int nu;
    int nv;
    
    // For convenience
    SEXP returnResult(SEXP d, SEXP u, SEXP v, SEXP nconv, SEXP niter)
    {
        return Rcpp::List::create(
            Rcpp::Named("d") = d,
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("nconv") = nconv,
            Rcpp::Named("niter") = niter
        );
    }
public:
    SVDsSym(int n_, int k_, int nu_, int nv_, int ncv_, MatOp *op_,
            double tol_ = 1e-10, int maxitr_ = 1000);
    Rcpp::List extract();
    ~SVDsSym() {}
};

#endif // SVDSSYM_H
