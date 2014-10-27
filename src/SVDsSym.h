#ifndef SVDSSYM_H
#define SVDSSYM_H

#include "EigsSym.h"

class SVDsSym: protected EigsSym
{
protected:
    int nu;
    int nv;
    
    // Converged singular values
    Rcpp::NumericVector dval;
    // Index of singular vectors associated with dval
    Rcpp::IntegerVector dind;
    
    // Extract eigenvectors of the matrix. In this case
    // they are also singular vectors
    Rcpp::RObject extractEigenvectors(int num);
    
    // For convenience
    Rcpp::List returnResult(SEXP d, SEXP u, SEXP v, SEXP nconv, SEXP niter)
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
    void compute();
    virtual Rcpp::List extract();
    ~SVDsSym() {}
};

#endif // SVDSSYM_H
