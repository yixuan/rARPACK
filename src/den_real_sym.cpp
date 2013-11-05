#include <RcppEigen.h>
#include "do_eigs.h"

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LDLT;
using Eigen::Lower;
using Eigen::Upper;
typedef Eigen::Map<Eigen::VectorXd> MapVec;


// Dense matrix-vector product
void den_sym_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    char *uplo = (char *) data;
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dsymv)(uplo, &n,
            &alpha, REAL(mat), &n,
            x_in, &one, &zero,
            y_out, &one);
}


// Data passed to den_sym*_mat_v_prod_shinv()
// L and U represent different solver using lower and upper triangles
typedef struct {
    const LDLT<MatrixXd> *solver;
    MapVec *x_vec;
    MapVec *y_vec;
} SymLData;

typedef struct {
    const LDLT<MatrixXd, Upper> *solver;
    MapVec *x_vec;
    MapVec *y_vec;
} SymUData;

// Shift-and-invert mode
void den_syml_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    SymLData *sldata = (SymLData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (sldata->x_vec) MapVec(x_in, n);
    new (sldata->y_vec) MapVec(y_out, n);
    *(sldata->y_vec) = (*(sldata->solver)).solve(*(sldata->x_vec));
}
void den_symu_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    SymUData *sudata = (SymUData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (sudata->x_vec) MapVec(x_in, n);
    new (sudata->y_vec) MapVec(y_out, n);
    *(sudata->y_vec) = (*(sudata->solver)).solve(*(sudata->x_vec));
}


RcppExport SEXP den_real_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r, SEXP sigma_scalar_r,
        SEXP workmode_scalar_r,
        SEXP lower_logical_r)
{
BEGIN_RCPP
    if(INTEGER(workmode_scalar_r)[0] == 1)
    {
        char uplo = LOGICAL(lower_logical_r)[0] ? 'L' : 'U';
        
        return do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                   which_string_r, ncv_scalar_r,
                   tol_scalar_r, maxitr_scalar_r,
                   retvec_logical_r, sigma_scalar_r,
                   workmode_scalar_r,
                   den_sym_mat_v_prod, &uplo);
    } else {
        int n = INTEGER(n_scalar_r)[0];
        
        // Map A_mat_r to Eigen matrix
        Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
        // Subtract the diagonal elements by sigma, i.e., A - sigma * I
        for(int i = 0; i < n; i++)
        {
            A(i, i) -= REAL(sigma_scalar_r)[0];
        }
        
        // Declare Eigen vectors (hey here I mean the C++ library Eigen,
        // not eigenvectors) that will be connected to workd in the iteration
        MapVec x_vec(NULL, n);
        MapVec y_vec(NULL, n);
        
        // If we use lower triangle of A
        if(LOGICAL(lower_logical_r)[0])
        {
            // Choleskey decomposition
            const LDLT<MatrixXd> cholA(A);
            // Data passed to den_syml_mat_v_prod_shinv()
            SymLData data;
            data.solver = &cholA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += REAL(sigma_scalar_r)[0];
            }
            
            return do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                       which_string_r, ncv_scalar_r,
                       tol_scalar_r, maxitr_scalar_r,
                       retvec_logical_r, sigma_scalar_r,
                       workmode_scalar_r,
                       den_syml_mat_v_prod_shinv, &data);
        } else {
            // Choleskey decomposition
            const LDLT<MatrixXd, Upper> cholA(A);
            // Data passed to den_symu_mat_v_prod_shinv()
            SymUData data;
            data.solver = &cholA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += REAL(sigma_scalar_r)[0];
            }

            return do_eigs_sym(A_mat_r, n_scalar_r, k_scalar_r,
                       which_string_r, ncv_scalar_r,
                       tol_scalar_r, maxitr_scalar_r,
                       retvec_logical_r, sigma_scalar_r,
                       workmode_scalar_r,
                       den_symu_mat_v_prod_shinv, &data);
        }
    }
    
    // Should not get here
    return R_NilValue;

END_RCPP
}
