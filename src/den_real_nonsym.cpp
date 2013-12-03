#include <RcppEigen.h>
#include "do_eigs.h"

using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::PartialPivLU;

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


// Dense matrix-vector product
// This function uses BLAS to calculte y_out = Mat * x_in
static void den_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                           int m, int n, void *data)
{
    char trans = 'N';
    double alpha = 1.0;
    int one = 1;
    double zero = 0.0;

    F77_CALL(dgemv)(&trans, &m, &n,
            &alpha, REAL(mat), &m,
            x_in, &one, &zero,
            y_out, &one);
}


// Data passed to denr_mat_v_prod_shinv()
typedef struct {
    const PartialPivLU<MatrixXd> *solver;
    MapVec *x_vec;
    MapVec *y_vec;
} DenRLUData; // dense, real, LU decomposition

// Data passed to denc_mat_v_prod_shinv()
typedef struct {
    const PartialPivLU<MatrixXcd> *solver;
    VectorXcd *x_vec;
    MapVec *y_vec;
} DenCLUData; // dense, complex, LU decomposition


// The two functions below support the Shift-and-invert mode of ARPACK
// They use LU decomposition to solve the linear equation
//                      Mat * y_out = x_in
// Equivalent to calculate y_out = Inv(Mat) * x_in
// denc_mat_v_prod_shinv() calculates y_out = Re(Inv(Mat)) * x_in
static void denr_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                                  int m, int n, void *data)
{
    DenRLUData *drludata = (DenRLUData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (drludata->x_vec) MapVec(x_in, n);
    new (drludata->y_vec) MapVec(y_out, m);
    (*(drludata->y_vec)).noalias() = (*(drludata->solver)).solve(*(drludata->x_vec));
}

static void denc_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                                  int m, int n, void *data)
{
    DenCLUData *dcludata = (DenCLUData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    (*(dcludata->x_vec)).real() = MapVec(x_in, n);
    new (dcludata->y_vec) MapVec(y_out, m);
    (*(dcludata->y_vec)).noalias() = (*(dcludata->solver)).solve(*(dcludata->x_vec)).real();
}



// Main function to solve dense, real, nonsymmetric eigen problems
RcppExport SEXP den_real_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                SEXP params_list_r,
                                SEXP sigmareal_logical_r)
{
BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    
    if(workmode == 1)
    {
        // If we don't need shift-and-invert mode
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,  
                              den_mat_v_prod, NULL);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // If sigma is real
            
            // Map A_mat_r to Eigen matrix
            MapMat A(REAL(A_mat_r), n, n);
        
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
        
            // Declare Eigen vectors
            // (hey here I mean the C++ library Eigen, not eigenvectors)
            // that will be connected to ARPACK
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denr_mat_v_prod_shinv()
            DenRLUData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  denr_mat_v_prod_shinv, &data);
        } else {
            // If sigma is complex
            
            // Map A_mat_r to Eigen matrix
            MapMat AR(REAL(A_mat_r), n, n);
            
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
        
            // Declare Eigen vectors that will be connected to ARPACK
            VectorXcd x_vec(n);
            x_vec.setZero();
            MapVec y_vec(NULL, n);
            
            // Data passed to denc_mat_v_prod_shinv()
            DenCLUData data;
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
} DenCLUData2;

// Shift-and-invert mode
static void denc_mat_v_prod_shinv2(SEXP mat, double *x_in, double *y_out,
                                   int m, int n, void *data)
{
    DenCLUData2 *dcludata = (DenCLUData2 *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (dcludata->x_vec) MapVec(x_in, n);
    new (dcludata->y_vec) MapVec(y_out, m);
    (*(dcludata->y_vec)).noalias() = (*(dcludata->solverA)).solve(*(dcludata->x_vec)) -
        (*(dcludata->solverS)).solve(*(dcludata->x_vec));
}

// Use real matrix operations to calculate Re(Inv(M))
// May be faster than den_real_nonsym() in standard optimization flags,
// but seems slower when -msse4 is used.
RcppExport SEXP den_real_nonsym2(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                 SEXP params_list_r,
                                 SEXP sigmareal_logical_r)
{
BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    
    if(workmode == 1)
    {
        // If we don't need shift-and-invert mode
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,  
                              den_mat_v_prod, NULL);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // If sigma is real
            
            // Map A_mat_r to Eigen matrix
            MapMat A(REAL(A_mat_r), n, n);
        
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
        
            // Declare Eigen vectors
            // (hey here I mean the C++ library Eigen, not eigenvectors)
            // that will be connected to ARPACK
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denr_mat_v_prod_shinv()
            DenRLUData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  denr_mat_v_prod_shinv, &data);
        } else {
            // If sigma is complex
            
            // Map A_mat_r to Eigen matrix
            MapMat A(REAL(A_mat_r), n, n);
        
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
        
            // Declare Eigen vectors that will be connected to ARPACK
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denc_mat_v_prod_shinv()
            DenCLUData2 data;
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
    
    if(workmode == 1)
    {
        // If we don't need shift-and-invert mode
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,  
                              den_mat_v_prod, NULL);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // If sigma is real
            
            // Map A_mat_r to Eigen matrix
            MapMat A(REAL(A_mat_r), n, n);
        
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
        
            // Declare Eigen vectors
            // (hey here I mean the C++ library Eigen, not eigenvectors)
            // that will be connected to ARPACK
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denr_mat_v_prod_shinv()
            DenRLUData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  denr_mat_v_prod_shinv, &data);
        } else {
            // If sigma is complex
            
            // Map A_mat_r to Eigen matrix
            MapMat A(REAL(A_mat_r), n, n);
        
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
        
            // Declare Eigen vectors that will be connected to ARPACK
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to denc_mat_v_prod_shinv()
            DenCLUData2 data;
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
