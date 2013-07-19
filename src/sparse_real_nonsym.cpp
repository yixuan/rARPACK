#include <RcppEigen.h>
#include "do_eigs.h"

typedef Eigen::MappedSparseMatrix<double> MapSpMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

// Data passed to matrix-vector product function
typedef struct {
    const MapSpMat *A;
    MapVec *x_vec;
    MapVec *y_vec;
} SparseData;

// Sparse matrix-vector product
void sparse_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                       int n, void *data)
{
    SparseData *spdata = (SparseData *) data;
    // first map x_in and y_out to x_vec and y_vec respectively,
    // and then do the matrix product y_vec <- A * x_vec
    new (spdata->x_vec) MapVec(x_in, n);
    new (spdata->y_vec) MapVec(y_out, n);
    (*(spdata->y_vec)).noalias() = *(spdata->A) * *(spdata->x_vec);
}


RcppExport SEXP sparse_real_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
        SEXP which_string_r, SEXP ncv_scalar_r,
        SEXP tol_scalar_r, SEXP maxitr_scalar_r,
        SEXP retvec_logical_r,
        SEXP sigmar_scalar_r, SEXP sigmai_scalar_r)
{
BEGIN_RCPP
    
    int n = INTEGER(n_scalar_r)[0];
    
    // map A_mat_r to Eigen sparse matrix
    const MapSpMat A(Rcpp::as<MapSpMat>(A_mat_r));
    // declare Eigen vectors (hey here I mean the C++ library Eigen,
    // not eigenvectors) that will be connected to workd in the iteration
    MapVec x_vec(NULL, n);
    MapVec y_vec(NULL, n);
    
    // Data passed to sparse_mat_v_prod()
    SparseData data;
    data.A = &A;
    data.x_vec = &x_vec;
    data.y_vec = &y_vec;
    
    return do_eigs_nonsym(A_mat_r, n_scalar_r, k_scalar_r,
                   which_string_r, ncv_scalar_r,
                   tol_scalar_r, maxitr_scalar_r,
                   retvec_logical_r,
                   sigmar_scalar_r, sigmai_scalar_r,
                   sparse_mat_v_prod, &data);

END_RCPP
}
