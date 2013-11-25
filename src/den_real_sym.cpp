#include <RcppEigen.h>
#include "do_eigs.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::PartialPivLU;

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


// Dense matrix-vector product
// This function uses BLAS to calculte y_out = Mat * x_in
// where Mat is a symmetric matrix
void den_sym_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                        int n, void *data)
{
    // *data is eigher 'L' or 'U'
    char *uplo = (char *) data;
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dsymv)(uplo, &n,
            &alpha, REAL(mat), &n,
            x_in, &one, &zero,
            y_out, &one);
}


// Data passed to den_sym_mat_v_prod_shinv()
typedef struct {
    const PartialPivLU<MatrixXd> *solver;
    MapVec *x_vec;
    MapVec *y_vec;
} SymLUData; // dense, real, symmetric, LU decomposition

// This functions supports the Shift-and-invert mode of ARPACK
// It uses LU decomposition to solve the linear equation
//                      Mat * y_out = x_in
// Equivalent to calculate y_out = Inv(Mat) * x_in
void den_sym_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                              int n, void *data)
{
    SymLUData *symludata = (SymLUData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (symludata->x_vec) MapVec(x_in, n);
    new (symludata->y_vec) MapVec(y_out, n);
    (*(symludata->y_vec)).noalias() = (*(symludata->solver)).solve(*(symludata->x_vec));
}



// Main function to solve dense, real, symmetric eigen problems
RcppExport SEXP den_real_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                             SEXP params_list_r,
                             SEXP lower_logical_r)
{
BEGIN_RCPP
    
    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    double sigma = as<double>(params_rcpp["sigma"]);
    
    if(workmode == 1)
    {
        // If we don't need shift-and-invert mode
        
        char uplo = LOGICAL(lower_logical_r)[0] ? 'L' : 'U';
        
        return do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                           params_list_r,
                           den_sym_mat_v_prod, &uplo);
    } else {
        
        int n = INTEGER(n_scalar_r)[0];
        
        // Map A_mat_r to Eigen matrix
        MapMat A(REAL(A_mat_r), n, n);
        
        // Subtract the diagonal elements by sigma, i.e., A - sigma * I
        for(int i = 0; i < n; i++)
        {
            A(i, i) -= sigma;
        }
        
        // LU decomposition
        const PartialPivLU<MatrixXd> luA(A);
        
        // Change A back
        for(int i = 0; i < n; i++)
        {
            A(i, i) += sigma;
        }
        
        // Declare Eigen vectors
        // (hey here I mean the C++ library Eigen, not eigenvectors)
        // that will be connected to ARPACK
        MapVec x_vec(NULL, n);
        MapVec y_vec(NULL, n);
        
        // Data passed to den_sym_mat_v_prod_shinv()
        SymLUData data;
        data.solver = &luA;
        data.x_vec = &x_vec;
        data.y_vec = &y_vec;
              
        return do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                           params_list_r,
                           den_sym_mat_v_prod_shinv, &data);
    }
    
    // Should not get here
    return R_NilValue;

END_RCPP
}
