#ifndef EIGSSYM_H
#define EIGSSYM_H

#include "Eigs.h"
#include "ARPACK.h"


// Implemented:
//     error()
//     warning()
//     aupd()
//     extract()
//
class EigsSym: public Eigs
{
protected:
    virtual void error(int stage, int errorcode);
    virtual void warning(int stage, int errorcode);

    // Store final results of eigenvectors
    // Conceptually it is an n * ncv matrix
    Rcpp::NumericMatrix eigV;
    // Store final results of eigenvalues
    // Conceptually it is a vector of length nev
    Rcpp::NumericVector eigd;

    virtual void aupd();
public:
    EigsSym(int n_, int nev_, int ncv_, MatOp *op_,
            const string & which_ = "LM", int workmode_ = 1, 
            char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual Rcpp::List extract(bool rvec = true);
    virtual ~EigsSym();
};


#endif // EIGSSYM_H

