#ifndef EIGS_H
#define EIGS_H

#include <Rcpp.h>
#include "MatOp.h"

using std::string;


// Need to be implemented:
//     error()
//     warning()
//     aupd()
//     extract()
//
class Eigs
{
protected:
    // Parameters to be set
    //
    // 'I' means standard eigenvalue problem,
    //     A * x = lambda * x
    // 'G' for generalized eigenvalue problem,
    //     A * x = lambda * B * x
    char bmat;
    // Dimension of A (n x n)
    int n;
    // Specify selection criteria
    // "LM": largest magnitude
    // "SM": smallest magnitude
    // "LR", "LI": largest real/imaginary part
    // "SR", "SI": smallest real/imaginary part
    // "LA": largest (algebraic) eigenvalues
    // "SA": smallest (algebraic) eigenvalues
    // "BE": half of each end
    string which;
    // Number of eigenvalues requested
    int nev;
    // Precision
    double tol;
    // Related to the algorithm, large ncv results in
    // faster convergence, but with greater memory use
    int ncv;
    // Maximum number of iterations
    int maxitr;
    // Workmode
    // workmode = 1: regular mode, we only need matrix product
    // workmode = 3: shift-and-invert mode,
    //               need to solve linear equation
    int workmode;

    // Matrix operation object
    MatOp *op;

    // Variables for computation in ARPACK
    //
    // Convergence flag
    int ido;
    // Warning and error flag in the first stage
    int info;
    // Warning and error flag in the second stage
    int ierr;
    // Control parameters
    int iparam[11];
    // Some pointers
    int ipntr[14];
    // Residual vector
    double *resid;
    // Working space and corresponding dimensions
    double *workd;
    int lworkl;
    double *workl;

    // stage = 1: _aupd
    // stage = 2: _eupd
    virtual void error(int stage, int errorcode) = 0;
    virtual void warning(int stage, int errorcode) = 0;
    
    // Generate initial residual vector
    void initResid();
    // Wrapper of _aupd
    virtual void aupd() = 0;
    // Check any error after update()
    void checkUpdateError();
    // Number of calls of _aupd(). Give an error if
    // extract() is called but updatecount == 0
    int updatecount;
public:
    // Constructor
    Eigs(int n_, int nev_, int ncv_, MatOp *op_,
         const string & which_ = "LM", int workmode_ = 1,
         char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    // _aupd step
    void update();
    // _eupd step, to extract results
    virtual Rcpp::List extract(bool rvec = true) = 0;
    // Destructor
    virtual ~Eigs();
};


#endif // EIGS_H

