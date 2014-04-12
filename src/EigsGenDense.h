#ifndef EIGSGENDENSE_H
#define EIGSGENDENSE_H

#include <RcppEigen.h>
#include "EigsGen.h"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::PartialPivLU;

typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;


// Implemented:
//     MultVector()
//     BindMatrix()
//
class EigsGenDense: public EigsGen
{
protected:
    // Pointer to the A matrix
    double *A_pntr;
    // Some constants for BLAS
    static const char BLAS_trans;
    static const double BLAS_alpha;
    static const int BLAS_one;
    static const double BLAS_zero;
    // Matrix inverse solver
    PartialPivLU<MatrixXd> solver;
    PartialPivLU<MatrixXcd> csolver;
    VectorXcd cx_vec;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;

public:
    EigsGenDense(int n_, int nev_, int ncv_,
                 const string & which_ = "LM", int workmode_ = 1,
                 double sigmar_ = 0, double sigmai_ = 0,
                 char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual void BindMatrix(SEXP mat_);
    virtual void MultVector(double *x_in, double *y_out);
    virtual void MultVectorShift(double *x_in, double *y_out);
    virtual ~EigsGenDense();
};

#endif // EIGSGENDENSE_H

