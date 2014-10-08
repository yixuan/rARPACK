#include "MatOp_function.h"

const char MatOp_function::BLAS_notrans = 'N';
const char MatOp_function::BLAS_trans = 'T';
const double MatOp_function::BLAS_alpha = 1.0;
const int MatOp_function::BLAS_one = 1;
const double MatOp_function::BLAS_zero = 0.0;

void MatOp_function::init()
{
    if(!canSolve)  return;

    if(m != n) return;

    // If sigma is real
    // if(fabs(sigmai) < 1e-17)
    // {
    //     // Subtract the diagonal elements by sigma, i.e., A - sigma * I
    //     for(int i = 0; i < n; i++)
    //     {
    //         A(i, i) -= sigmar;
    //     }
    //     // LU decomposition
    //     solver.compute(A);
    //     // Change A back
    //     for(int i = 0; i < n; i++)
    //     {
    //         A(i, i) += sigmar;
    //     }
    // } else {
    //     MatrixXcd cA = A.cast< std::complex<double> >();
    //     // Subtract the diagonal elements by sigma, i.e., A - sigma * I
    //     for(int i = 0; i < n; i++)
    //     {
    //         cA(i, i) -= std::complex<double>(sigmar, sigmai);
    //     }
    //     // LU decomposition
    //     csolver.compute(cA);

    //     cx_vec.resize(n);
    //     cx_vec.setZero();
    // }
}

MatOp_function::MatOp_function(Rcpp::Function FUN_function_r, int m_, int n_,
                 bool needSolve_):
    MatOp(m_, n_, 0, 0, true, needSolve_),
    FUN_function_(FUN_function_r)
{
    init();
}

void MatOp_function::prod(double *x_in, double *y_out)
{
    Rcpp::NumericVector y(n);
    Rcpp::NumericVector x(n, 1.0);
    for(int i = 0; i<n; i++){
        x[i] = *(x_in + i);
    }
    y = FUN_function_(x);
    for(int i = 0; i<n; i++){
        *(y_out + i) = y[i];
    }
}

void MatOp_function::tprod(double *x_in, double *y_out)
{
    Rcpp::NumericVector y(n);
    Rcpp::NumericVector x(n, 1.0);
    for(int i = 0; i<n; i++){
        x[i] = *(x_in + i);
    }
    y = FUN_function_(x);
    for(int i = 0; i<n; i++){
        *(y_out + i) = y[i];
    }
}

// void MatOp_function::shiftSolve(double *x_in, double *y_out)
// {
//     if(m != n)
//         Rcpp::stop("matrix must be square");
//     if(!canSolve)
//         Rcpp::stop("this matrix doesn't support solving linear equation");

//     if(fabs(sigmai) < 1e-17)
//     {
//         new (&x_vec) MapVec(x_in, n);
//         new (&y_vec) MapVec(y_out, n);
//         y_vec = solver.solve(x_vec);
//     } else {
//         cx_vec.real() = MapVec(x_in, n);
//         new (&y_vec) MapVec(y_out, n);
//         y_vec = csolver.solve(cx_vec).real();
//     }
// }
