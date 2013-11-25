#include <RcppEigen.h>
#include "do_eigs.h"

using Rcpp::as;
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
                                SEXP params_list_r,
                                SEXP sigmareal_logical_r)
{
BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    
    // If we don't need shift-and-invert mode
    if(workmode == 1)
    {
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,  
                              den_mat_v_prod, NULL);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= sigmar;
            }
            
            // LU decomposition
            const PartialPivLU<MatrixXd> luA(A);
        
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += sigmar;
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
                                  params_list_r,
                                  denr_mat_v_prod_shinv, &data);
        } else { 
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> AR(REAL(A_mat_r), n, n);
            
            //MatrixXcd A(n, n);
            //A.real() = AR;
            //A.imag().setZero();
            MatrixXcd A = AR.cast< std::complex<double> >();
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= std::complex<double>(sigmar, sigmai);
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
                                  params_list_r,
                                  denc_mat_v_prod_shinv, &data);
        }
    }
    
    // Should not get here
    return R_NilValue;
    

END_RCPP
}





// Alternative methods
//
//
//
//
typedef struct {
    const PartialPivLU<MatrixXd> *solverA;
    const PartialPivLU<MatrixXd> *solverS;
    MapVec *x_vec;
    MapVec *y_vec;
} DenCData2;

// Shift-and-invert mode
void denc_mat_v_prod_shinv2(SEXP mat, double *x_in, double *y_out,
                            int n, void *data)
{
    DenCData2 *dcdata = (DenCData2 *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (dcdata->x_vec) MapVec(x_in, n);
    new (dcdata->y_vec) MapVec(y_out, n);
    *(dcdata->y_vec) = (*(dcdata->solverA)).solve(*(dcdata->x_vec)) -
        (*(dcdata->solverS)).solve(*(dcdata->x_vec));
}

// Using real matrix operations to calculate Re(Inv(M))
// May be faster than den_real_nonsym() in standard optimization flags,
// but seems slower when -msse4 is used.
RcppExport SEXP den_real_nonsym2(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                 SEXP params_list_r,
                                 SEXP sigmareal_logical_r)
{
BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    
    // If we don't need shift-and-invert mode
    if(workmode == 1)
    {
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,  
                              den_mat_v_prod, NULL);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= sigmar;
            }
            
            // LU decomposition
            const PartialPivLU<MatrixXd> luA(A);
        
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += sigmar;
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
                                  params_list_r,
                                  denr_mat_v_prod_shinv, &data);
        } else {
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= sigmar;
            }
            
            // sigma = a + b * i
            // Aa = A - a * I
            // A - sigma * I = Aa - b * I * i
            // inv(A - sigma * I) = inv(Aa) - inv(b^(-2) * Aa^3 + Aa)
            //                    = inv(Aa) - inv(S)
            // S <- A / b^2
            // S <- S * A + I
            // S <- S * A
            double b2 = sigmai * sigmai;
            MatrixXd S(A);

            // Use Eigen
            S /= b2;
            S = S * A;
            for(int i = 0; i < n; i++)
            {
                S(i, i) += 1;
            }
            S = S * A;
            
            // LU decomposition
            const PartialPivLU<MatrixXd> luA(A);
            const PartialPivLU<MatrixXd> luS(S);
        
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += sigmar;
            }
        
            // Declare Eigen vectors (hey here I mean the C++ library Eigen,
            // not eigenvectors) that will be connected to workd in the iteration
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denc_mat_v_prod_shinv()
            DenCData2 data;
            data.solverA = &luA;
            data.solverS = &luS;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  denc_mat_v_prod_shinv2, &data);
        }
    }
    
    // Should not get here
    return R_NilValue;
    

END_RCPP
}

// Optimized using BLAS. Will be fast using a highly optimized BLAS.
RcppExport SEXP den_real_nonsym3(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                 SEXP params_list_r,
                                 SEXP sigmareal_logical_r)
{
BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    
    // If we don't need shift-and-invert mode
    if(workmode == 1)
    {
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,  
                              den_mat_v_prod, NULL);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= sigmar;
            }
            
            // LU decomposition
            const PartialPivLU<MatrixXd> luA(A);
        
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += sigmar;
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
                                  params_list_r,
                                  denr_mat_v_prod_shinv, &data);
        } else {
            // Map A_mat_r to Eigen matrix
            Map<MatrixXd> A(REAL(A_mat_r), n, n);
        
            // Subtract the diagonal elements by sigma, i.e., A - sigma * I
            for(int i = 0; i < n; i++)
            {
                A(i, i) -= sigmar;
            }
            
            // sigma = a + b * i
            // Aa = A - a * I
            // A - sigma * I = Aa - b * I * i
            // inv(A - sigma * I) = inv(Aa) - inv(b^(-2) * Aa^3 + Aa)
            //                    = inv(Aa) - inv(S)
            // S <- A / b^2
            // S <- S * A + I
            // S <- S * A
            double b2 = sigmai * sigmai;
            MatrixXd S(A);
            // Use BLAS
            // C <- alpha * A * B + beta * C
            // F77_CALL(dgemm)(&trans, &trans, &n, &n, &n,
            //     &alpha, A, &n, B, &n,
            //     &beta, C, &n);
            
            char trans = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            MatrixXd tmp(n, n);
            // tmp <- A * A
            F77_CALL(dgemm)(&trans, &trans, &n, &n, &n,
                &alpha, &A(0, 0), &n, &A(0, 0), &n,
                &beta, &tmp(0, 0), &n);
            // S <- alpha * tmp * A + S
            alpha = 1.0 / b2;
            beta = 1.0;
            F77_CALL(dgemm)(&trans, &trans, &n, &n, &n,
                &alpha, &tmp(0, 0), &n, &A(0, 0), &n,
                &beta, &S(0, 0), &n);
            
            // LU decomposition
            const PartialPivLU<MatrixXd> luA(A);
            const PartialPivLU<MatrixXd> luS(S);
        
            // Change A back
            for(int i = 0; i < n; i++)
            {
                A(i, i) += sigmar;
            }
        
            // Declare Eigen vectors (hey here I mean the C++ library Eigen,
            // not eigenvectors) that will be connected to workd in the iteration
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denc_mat_v_prod_shinv()
            DenCData2 data;
            data.solverA = &luA;
            data.solverS = &luS;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  denc_mat_v_prod_shinv2, &data);
        }
    }
    
    // Should not get here
    return R_NilValue;
    

END_RCPP
}