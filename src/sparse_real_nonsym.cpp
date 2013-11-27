#include <RcppEigen.h>
#include "do_eigs.h"

using Rcpp::as;
using Eigen::VectorXi;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
using Eigen::SparseLU;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix< std::complex<double> > SpCMat;
typedef Eigen::MappedSparseMatrix<double> MapSpMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

// Data passed to sparse_mat_v_prod()
typedef struct {
    const MapSpMat *A;
    MapVec *x_vec;
    MapVec *y_vec;
} SparseData;

// Data passed to sparser_mat_v_prod_shinv()
typedef struct {
    SparseLU<SpMat> *solver;
    MapVec *x_vec;
    MapVec *y_vec;
} SparseRLUData; // sparse, real, LU decomposition

// Data passed to sparsec_mat_v_prod_shinv()
typedef struct {
    SparseLU<SpCMat> *solver;
    VectorXcd *x_vec;
    MapVec *y_vec;
} SparseCLUData; // sparse, complex, LU decomposition


// Sparse matrix-vector product
// This function calcultes y_out = Mat * x_in
// where Mat is a sparse matrix
void sparse_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                       int m, int n, void *data)
{
    SparseData *spdata = (SparseData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then do the matrix product y_vec <- A * x_vec
    new (spdata->x_vec) MapVec(x_in, n);
    new (spdata->y_vec) MapVec(y_out, m);
    (*(spdata->y_vec)).noalias() = *(spdata->A) * *(spdata->x_vec);
}

// The two functions below support the Shift-and-invert mode of ARPACK
// They use sparse LU decomposition to solve the linear equation
//                      Mat * y_out = x_in
// where Mat is a sparse matrix
// Equivalent to calculate y_out = Inv(Mat) * x_in
// sparsec_mat_v_prod_shinv() calculates y_out = Re(Inv(Mat)) * x_in
void sparser_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                              int m, int n, void *data)
{  
    SparseRLUData *sprludata = (SparseRLUData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    new (sprludata->x_vec) MapVec(x_in, n);
    new (sprludata->y_vec) MapVec(y_out, m);
    (*(sprludata->y_vec)).noalias() = (*(sprludata->solver)).solve(*(sprludata->x_vec));
}

void sparsec_mat_v_prod_shinv(SEXP mat, double *x_in, double *y_out,
                              int m, int n, void *data)
{  
    SparseCLUData *spcludata = (SparseCLUData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then solve the linear equation Mat * y_out = x_in
    (*(spcludata->x_vec)).real() = MapVec(x_in, n);
    new (spcludata->y_vec) MapVec(y_out, m);
    (*(spcludata->y_vec)).noalias() = (*(spcludata->solver)).solve(*(spcludata->x_vec)).real();
}



// Main function to solve sparse, real, nonsymmetric eigen problems
RcppExport SEXP sparse_real_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                   SEXP params_list_r,
                                   SEXP sigmareal_logical_r)
{
BEGIN_RCPP
    
    Rcpp::List params_rcpp(params_list_r);
    int workmode = as<int>(params_rcpp["workmode"]);
    int n = INTEGER(n_scalar_r)[0];
    
    if(workmode == 1)
    {
        // If we don't need shift-and-invert mode
        
        // Map A_mat_r to Eigen sparse matrix
        const MapSpMat A(as<MapSpMat>(A_mat_r));
  
        // Declare Eigen vectors
        // (hey here I mean the C++ library Eigen, not eigenvectors)
        // that will be connected to ARPACK
        MapVec x_vec(NULL, n);
        MapVec y_vec(NULL, n);
    
        // Data passed to sparse_mat_v_prod()
        SparseData data;
        data.A = &A;
        data.x_vec = &x_vec;
        data.y_vec = &y_vec;
    
        return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                              params_list_r,
                              sparse_mat_v_prod, &data);
    } else {
        
        double sigmar = as<double>(params_rcpp["sigmar"]);
        double sigmai = as<double>(params_rcpp["sigmai"]);
        int n = INTEGER(n_scalar_r)[0];
        
        if(LOGICAL(sigmareal_logical_r)[0])
        {
            // If sigma is real
            
            // Map A_mat_r to Eigen sparse matrix
            const MapSpMat A(as<MapSpMat>(A_mat_r));
        
            // Create a sparse identity matrix
            SpMat I(n, n);
            I.setIdentity();
            
            // Sparse LU decomposition
            SparseLU<SpMat> luA;
            luA.compute(A - sigmar * I);
            if(luA.info() != Eigen::Success) {
                ::Rf_error("Eigen: LU decomposition of A - sigma * I failed");
            }

            // Declare Eigen vectors that will be connected to ARPACK
            MapVec x_vec(NULL, n);
            MapVec y_vec(NULL, n);
            
            // Data passed to sparser_mat_v_prod_shinv()
            SparseRLUData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  sparser_mat_v_prod_shinv, &data);
        } else {
            // If sigma is complex
            
            // Map A_mat_r to Eigen matrix
            const MapSpMat AR(as<MapSpMat>(A_mat_r));
            
            SpCMat A = AR.cast< std::complex<double> >();
            
            // Create a sparse identity matrix (1 + 0i on diagonal)
            SpCMat I(n, n);
            I.setIdentity();
        
            // Sparse LU decomposition
            SparseLU<SpCMat> luA;
            luA.compute(A - std::complex<double>(sigmar, sigmai) * I);
            if(luA.info() != Eigen::Success) {
                ::Rf_error("Eigen: LU decomposition of A - sigma * I failed");
            }
        
            // Declare Eigen vectors that will be connected to ARPACK
            VectorXcd x_vec(n);
            x_vec.setZero();
            MapVec y_vec(NULL, n);
            
            // Data passed to sparsec_mat_v_prod_shinv()
            SparseCLUData data;
            data.solver = &luA;
            data.x_vec = &x_vec;
            data.y_vec = &y_vec;
            
            return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                                  params_list_r,
                                  sparsec_mat_v_prod_shinv, &data);
        }
    }
    
    // Should not get here
    return R_NilValue;


END_RCPP
}
