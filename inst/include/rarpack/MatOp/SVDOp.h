#ifndef SVDOP_H
#define SVDOP_H

#include <algorithm>
#include "MatProd.h"

class SVDTallOp: public MatProd
{
private:
    MapProd *op;
    const int dim;
    double *work;
public:
    SVDTallOp(MatProd *op_) :
        op(op_),
        dim(std::min(op->rows(), op->cols())),
        work(new double[op->rows()])
    {}

    ~SVDTallOp()
    {
        delete [] work;
    }

    int rows() { return dim; }
    int cols() { return dim; }

    // y_out = A'A * x_in
    void perform_op(double *x_in, double *y_out)
    {
        op->perform_op(x_in, tmp);
        op->perform_tprod(tmp, y_out);
    }
};

class SVDWideOp: public MatProd
{
private:
    MapProd *op;
    const int dim;
    double *work;
public:
    SVDWideOp(MatProd *op_) :
        op(op_),
        dim(std::min(op->rows(), op->cols())),
        work(new double[op->cols()])
    {}

    ~SVDWideOp()
    {
        delete [] work;
    }

    int rows() { return dim; }
    int cols() { return dim; }

    // y_out = AA' * x_in
    void perform_op(double *x_in, double *y_out)
    {
        op->perform_tprod(x_in, tmp);
        op->perform_op(tmp, y_out);
    }
};

#endif // SVDOP_H
