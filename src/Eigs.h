#ifndef EIGS_H
#define EIGS_H

#include <Rcpp.h>
using std::string;


// Need to be implemented:
//     Error()
//     Warning()
//     MultVector()
//     MultVectorShift()
//     BindMatrix()
//     AllocMem()
//     Update()
//     Extract()
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
    int workmode;
    // Shift parameter, real part
    double sigmar;
    // Shift parameter, imaginary part
    double sigmai;


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

    // Flag to indicate whether matrix has been linked
    bool matrix_linked;

    // stage = 1: _aupd
    // stage = 2: _eupd
    virtual void Error(int stage, int errorcode) = 0;
    virtual void Warning(int stage, int errorcode) = 0;
    
    // Generate initial residual vector
    void InitResid();

    // Whether the matrix has been linked to internal structure
    bool MatrixLinked() { return matrix_linked; }
    bool MatrixLinked(bool linked) { matrix_linked = linked;
                                     return matrix_linked; }
public:
    // Constructor
    Eigs(int n_, int nev_, int ncv_,
         const string & which_ = "LM", int workmode_ = 1,
         double sigmar_ = 0, double sigmai_ = 0,
         char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    // Map matrix in R to internal structures
    virtual void BindMatrix(SEXP mat_) = 0;
    // y_out = A * x_in
    virtual void MultVector(double *x_in, double *y_out) = 0;
    // y_out = inv(A - sigma * I) * x_in
    virtual void MultVectorShift(double *x_in, double *y_out) = 0;
    // _aupd step
    virtual void Update() = 0;
    // _eupd step, to extract results
    virtual Rcpp::List Extract(bool rvec = true) = 0;
    // Destructor
    virtual ~Eigs();
};


#endif // EIGS_H

