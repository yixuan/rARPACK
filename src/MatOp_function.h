#ifndef MATOP_FUNCTION_H
#define MATOP_FUNCTION_H

#include <Rcpp.h>
#include "MatOp.h"

using Rcpp::as;

// Operations on "function" class
class MatOp_function : public MatOp
{
private:
    // The R function applied to a vector
    Rcpp::Function FUN;
    // Arguments passed to FUN
    Rcpp::List args;
public:
    // Constructor
    MatOp_function(Rcpp::Function FUN_, Rcpp::List args_, int n_);
    // y_out = FUN(x_in)
    void prod(double *x_in, double *y_out);
    // Destructor
    virtual ~MatOp_function() {}
};


#endif // MATOP_FUNCTION_H
