#ifndef EIGSGEN_H
#define EIGSGEN_H

#include <RcppEigen.h>
#include "Eigs.h"
#include "ARPACK.h"


// Implemented:
//     error()
//     warning()
//     aupd()
//     eupd()
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
    double *V;
    // Store final results of eigenvalues,
    // both real part and imaginary part
    // Conceptually they are vectors of length nev + 1
    double *dr;
    double *di;

    virtual void aupd();
    virtual void eupd();

    // eigenvalue with positive imaginary part of a 2 by 2 matrix
    // [a  b]
    // [c  d]
    // when we know it has complex eigenvalues
    static std::complex<double> eigenvalue2x2(const double &a,
        const double &b, const double &c, const double &d);
    static void eigenvalueSchur(const Eigen::MatrixXd &Rm,
                                Eigen::VectorXcd &result);
    static void findMatchedIndex(const Eigen::VectorXcd &target,
                                 const Eigen::VectorXcd &collection,
                                 Eigen::VectorXi &result);
    void recomputeH();
    void transformEigenvalues(Eigen::VectorXcd &evals);
    
    static void sortDesc(Eigen::VectorXcd &values);
    static Eigen::VectorXi sortDescWithOrder(Eigen::VectorXcd &values);
public:
    EigsGen(int n_, int nev_, int ncv_, MatOp *op_,
            const std::string & which_ = "LM", int workmode_ = 1, 
            char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual Rcpp::List extract();
    virtual ~EigsGen();
};


#endif // EIGSGEN_H

