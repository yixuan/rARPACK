#ifndef EIGSGEN_H
#define EIGSGEN_H

#include "Eigs.h"
#include "do_eigs.h"


// Implemented:
//     Error()
//     Warning()
//     Update()
//     Extract()
//
// Need to be implemented:
//     MultVector()
//     MultVectorShift()
//     BindMatrix()
//
class EigsGen: public Eigs
{
protected:
    virtual void Error(int stage, int errorcode);
    virtual void Warning(int stage, int errorcode);
    // Check any error after the _aupd step
    void CheckUpdateError();
    int updatecount;
    // Working space
    double *workv;

    // Store final results of eigenvectors
    // Conceptually it is an n * ncv matrix
    Rcpp::NumericMatrix eigV;
    // Store final results of eigenvalues,
    // both real part and imaginary part
    // Conceptually they are vectors of length nev + 1
    Rcpp::NumericVector eigdr;
    Rcpp::NumericVector eigdi;
public:
    EigsGen(int n_, int nev_, int ncv_,
            const string & which_ = "LM", int workmode_ = 1, 
            double sigmar_ = 0, double sigmai_ = 0,
            char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual void Update();
    virtual Rcpp::List Extract(bool rvec = true);
    virtual ~EigsGen();
};


#endif // EIGSGEN_H

