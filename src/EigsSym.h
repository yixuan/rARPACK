#ifndef EIGSSYM_H
#define EIGSSYM_H

#include "Eigs.h"
#include "do_eigs.h"


// Implemented:
//     Error()
//     Warning()
//     AllocMem()
//     Update()
//     Extract()
//
// Need to be implemented:
//     MultVector()
//     MultVectorShift()
//     BindMatrix()
//
class EigsSym: public Eigs
{
protected:
    virtual void Error(int stage, int errorcode);
    virtual void Warning(int stage, int errorcode);
    // Check any error after the _aupd step
    void CheckUpdateError();
    int updatecount;

    // Store final results of eigenvectors
    // Conceptually it is an n * ncv matrix
    Rcpp::NumericMatrix eigV;
    // Store final results of eigenvalues
    // Conceptually it is a vector of length nev
    Rcpp::NumericVector eigd;
public:
    EigsSym(int n_, int nev_, int ncv_,
            const string & which_ = "LM", int workmode_ = 1, 
            double sigmar_ = 0,
            char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    void AllocMem();
    virtual void Update();
    virtual Rcpp::List Extract(bool rvec = true);
    virtual ~EigsSym();
};


#endif // EIGSSYM_H

