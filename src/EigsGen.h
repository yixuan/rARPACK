#ifndef EIGSGEN_H
#define EIGSGEN_H

#include "Eigs.h"
#include "ARPACK.h"


// Implemented:
//     error()
//     warning()
//     aupd()
//     extract()
//
class EigsGen: public Eigs
{
protected:
    virtual void error(int stage, int errorcode);
    virtual void warning(int stage, int errorcode);

    // Working space, unique to eigs() with general matrices
    double *workv;

    // Store final results of eigenvectors
    // Conceptually it is an n * ncv matrix
    Rcpp::NumericMatrix eigV;
    // Store final results of eigenvalues,
    // both real part and imaginary part
    // Conceptually they are vectors of length nev + 1
    Rcpp::NumericVector eigdr;
    Rcpp::NumericVector eigdi;

    virtual void aupd();
public:
    EigsGen(int n_, int nev_, int ncv_, MatOp *op_,
            const std::string & which_ = "LM", int workmode_ = 1, 
            char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual Rcpp::List extract(bool rvec = true);
    virtual ~EigsGen();
};


#endif // EIGSGEN_H

