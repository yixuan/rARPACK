#include "MatOp_function.h"

MatOp_function::MatOp_function(Rcpp::Function FUN_, int n_,
                               Rcpp::List args_) :
    MatOp(n_, n_, 0, 0, false, false),
    FUN(FUN_), args(args_)
{
}

void MatOp_function::prod(double *x_in, double *y_out)
{
    Rcpp::NumericVector x(n);
    std::copy(x_in, x_in + n, x.begin());

    Rcpp::NumericVector y = FUN(x, args);
    if(y.length() != n)
        Rcpp::stop("the provided function should return n elements");
    
    std::copy(y.begin(), y.end(), y_out);
}
