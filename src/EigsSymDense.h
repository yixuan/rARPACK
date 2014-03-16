#ifndef EIGSSYMDENSE_H
#define EIGSSYMDENSE_H

#include <RcppEigen.h>
#include "EigsSym.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::PartialPivLU;

typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;


// Implemented:
//     MultVector()
//     BindMatrix()
//
class EigsSymDense: public EigsSym
{
protected:
    // Pointer to the A matrix
    double *A_pntr;
    // Which part of A should be used
    // 'L' or 'U'
    char uplo;
    // Some constants for BLAS
    static const double BLAS_alpha;
    static const int BLAS_one;
    static const double BLAS_zero;
    // Matrix inverse solver
    PartialPivLU<MatrixXd> solver;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;

public:
    EigsSymDense(int n_, int nev_, int ncv_,
                 const string & which_ = "LM", int workmode_ = 1,
                 double sigmar_ = 0,
                 char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual void BindMatrix(SEXP mat_);
    virtual void BindMatrix(SEXP mat_, char uplo_);
    virtual void MultVector(double *x_in, double *y_out);
    virtual void MultVectorShift(double *x_in, double *y_out);
    virtual ~EigsSymDense();
};

#endif // EIGSSYMDENSE_H

