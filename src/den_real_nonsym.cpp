#include <RcppEigen.h>
#include "do_eigs.h"

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::PartialPivLU;
typedef Eigen::Map<Eigen::VectorXd> MapVec;


// Dense matrix-vector product
void den_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    char trans = 'N';
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dgemv)(&trans, &n, &n,
            &alpha, REAL(mat), &n,
            x_in, &one, &zero,
            y_out, &one);
}


// Data passed to denr_mat_v_prod_shinv()
typedef struct {
    const PartialPivLU<MatrixXd> *solver;
    MapVec *x_vec;
    MapVec *y_vec;
} DenRData;

typedef struct {
    const PartialPivLU<MatrixXcd> *solver;
    VectorXcd *x_vec;
    MapVec *y_vec;
} DenCData;


// Shift-and-invert mode
void denr_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    DenRData *drdata = (DenRData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (drdata->x_vec) MapVec(x_in, n);
    new (drdata->y_vec) MapVec(y_out, n);
    *(drdata->y_vec) = (*(drdata->solver)).solve(*(drdata->x_vec));
}

void denc_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                    int n, void *data)
{
    DenCData *dcdata = (DenCData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    (*(dcdata->x_vec)).real() = MapVec(x_in, n);
    new (dcdata->y_vec) MapVec(y_out, n);
    *(dcdata->y_vec) = (*(dcdata->solver)).solve(*(dcdata->x_vec)).real();
}



RcppExport SEXP den_real_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r,
        SEXP sigmar_scalar_r, SEXP sigmai_scalar_r,
        SEXP workmode_scalar_r, SEXP sigmareal_logical_r)
{
BEGIN_RCPP
    // If we don't need shift-and-invert mode
    if(INTEGER(workmode_scalar_r)[0] == 1)
    {
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                   which_string_r, ncv_scalar_r,
                   tol_scalar_r, maxitr_scalar_r,
                   retvec_logical_r,
                   sigmar_scalar_r, sigmai_scalar_r,
                   workmode_scalar_r,
                   den_mat_v_prod, NULL);
    } else {
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            int n = INTEGER(n_scalar_r)[0];
        
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= REAL(sigmar_scalar_r)[0];
            }
            
            // LU decomposition
            const PartialPivLU<MatrixXd> luA(A);
        
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += REAL(sigmar_scalar_r)[0];
            }
        
            // Declare Eigen vectors (hey here I mean the C++ library Eigen,
            // not eigenvectors) that will be connected to workd in the iteration
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denr_mat_v_prod_shinv()
            DenRData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                       which_string_r, ncv_scalar_r,
                       tol_scalar_r, maxitr_scalar_r,
                       retvec_logical_r,
                       sigmar_scalar_r, sigmai_scalar_r,
                       workmode_scalar_r,
                       denr_mat_v_prod_shinv, &data);
        } else {
            int n = INTEGER(n_scalar_r)[0];
        
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> AR(REAL(A_mat_r), n, n);
            
            MatrixXcd A(n, n);
            A.real() = AR;
            A.imag().setZero();
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= std::complex<double>(REAL(sigmar_scalar_r)[0],
                                                REAL(sigmai_scalar_r)[0]);
            }
            
            // LU decomposition
            const PartialPivLU<MatrixXcd> luA(A);
        
            // Declare Eigen vectors (hey here I mean the C++ library Eigen,
            // not eigenvectors) that will be connected to workd in the iteration
            VectorXcd x_vec(n);
            x_vec.imag().setZero();
            MapVec y_vec(NULL, n);
            
            // Data passed to denc_mat_v_prod_shinv()
            DenCData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                       which_string_r, ncv_scalar_r,
                       tol_scalar_r, maxitr_scalar_r,
                       retvec_logical_r,
                       sigmar_scalar_r, sigmai_scalar_r,
                       workmode_scalar_r,
                       denc_mat_v_prod_shinv, &data);
        }
    }
    
    // Should not get here
    return R_NilValue;
    

END_RCPP
}

