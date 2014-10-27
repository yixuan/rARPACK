#ifndef SVDSGEN_H
#define SVDSGEN_H

#include "SVDsSym.h"

class SVDsGenTall: public SVDsSym
{
protected:
    int nrow;
    int ncol;
    // Workspace, used to calculate OP(x) = A' * A * x
    double *workspace;
    
    // OP(x) = A' * A * x
    virtual void matOp(double *x_in, double *y_out)
    {
        op->prod(x_in, workspace);
        op->tprod(workspace, y_out);
    }
    // matProd = op->prod
    virtual void matProd(double *x_in, double *y_out)
    {
        op->prod(x_in, y_out);
    }
    
    // Compute U, the left singular vectors
    Rcpp::RObject computeU(int num);
public:
    SVDsGenTall(int nrow_, int ncol_, int k_, int nu_, int nv_,
                int ncv_, MatOp *op_,
                double tol_ = 1e-10, int maxitr_ = 1000) :
        SVDsSym(ncol_, k_, nu_, nv_, ncv_, op_, tol_, maxitr_),
        nrow(nrow_), ncol(ncol_)
    {
        workspace = new double[nrow_];
    }
    virtual Rcpp::List extract();
    virtual ~SVDsGenTall() { delete [] workspace; }
};


class SVDsGenWide: public SVDsGenTall
{
private:    
    // OP(x) = A * A' * x
    void matOp(double *x_in, double *y_out)
    {
        op->tprod(x_in, workspace);
        op->prod(workspace, y_out);
    }
    // matProd = op->tprod
    void matProd(double *x_in, double *y_out)
    {
        op->tprod(x_in, y_out);
    }
    
    // Compute V, the right singular vectors
    Rcpp::RObject computeV(int num) { return computeU(num); }
public:
    SVDsGenWide(int nrow_, int ncol_, int k_, int nu_, int nv_,
                int ncv_, MatOp *op_,
                double tol_ = 1e-10, int maxitr_ = 1000) :
        SVDsGenTall(ncol_, nrow_, k_, nu_, nv_, ncv_, op_, tol_, maxitr_)
    {
        nrow = nrow_;
        ncol = ncol_;
    }
    Rcpp::List extract();
    ~SVDsGenWide() {}
};

#endif // SVDSGEN_H
