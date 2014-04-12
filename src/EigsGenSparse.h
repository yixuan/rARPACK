#ifndef EIGSGENSPARSE_H
#define EIGSGENSPARSE_H

#include <RcppEigen.h>
#include "EigsGen.h"

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::SparseLU;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix< std::complex<double> > SpCMat;
typedef Eigen::MappedSparseMatrix<double> MapSpMat;
typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;


// Implemented:
//     MultVector()
//     BindMatrix()
//
class EigsGenSparse: public EigsGen
{
protected:
    // Sparse matrix structure
    MapSpMat A;
    // Matrix inverse solver
    SparseLU<SpMat> solver;
    SparseLU<SpCMat> csolver;
    VectorXcd cx_vec;
    // Mapped vector
    MapVec x_vec;
    MapVec y_vec;

public:
    EigsGenSparse(int n_, int nev_, int ncv_,
                 const string & which_ = "LM", int workmode_ = 1,
                 double sigmar_ = 0, double sigmai_ = 0,
                 char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual void BindMatrix(SEXP mat_);
    virtual void MultVector(double *x_in, double *y_out);
    virtual void MultVectorShift(double *x_in, double *y_out);
    virtual ~EigsGenSparse();
};

#endif // EIGSGENSPARSE_H

