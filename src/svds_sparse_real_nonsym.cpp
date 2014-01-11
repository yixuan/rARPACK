#include <RcppEigen.h>
#include "do_eigs.h"

using Rcpp::as;

typedef Eigen::MappedSparseMatrix<double> MapSpMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;


// Data passed to sparse_mat_v_prod() and sparse_mat_t_v_prod()
typedef struct {
    const MapSpMat *A;
    MapVec *x_vec;
    MapVec *y_vec;
} SparseData;

// Calculate y_out = Mat * x_in for sparse matrix
static void sparse_mat_v_prod(SEXP mat, double *x_in, double *y_out,
                              int m, int n, void *data)
{
    SparseData *spdata = (SparseData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then do the matrix product y_vec <- Mat * x_vec
    new (spdata->x_vec) MapVec(x_in, n);
    new (spdata->y_vec) MapVec(y_out, m);
    (*(spdata->y_vec)).noalias() = *(spdata->A) * *(spdata->x_vec);
}

// Calculate y_out = Mat' * x_in for sparse matrix
static void sparse_mat_t_v_prod(SEXP mat, double *x_in, double *y_out,
                                int m, int n, void *data)
{
    SparseData *spdata = (SparseData *) data;
    // First map x_in and y_out to x_vec and y_vec respectively,
    // and then do the matrix product y_vec <- Mat' * x_vec
    new (spdata->x_vec) MapVec(x_in, m);
    new (spdata->y_vec) MapVec(y_out, n);
    (*(spdata->y_vec)).noalias() = (*(spdata->A)).transpose() * *(spdata->x_vec);
}


// Main function to calculate truncated SVD for
// sparse, real, nonsymmetric matrices
RcppExport SEXP sparse_real_nonsym_svd(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                                       SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                                       SEXP params_list_r)
{
BEGIN_RCPP

    // Map A_mat_r to Eigen sparse matrix
    const MapSpMat A(as<MapSpMat>(A_mat_r));
  
    // Declare in/out vectors
    // Initial pointer and length are not important
    MapVec x_vec(NULL, 1);
    MapVec y_vec(NULL, 1);
    
    // Data passed to sparse_mat_v_prod() and sparse_mat_t_v_prod()
    SparseData data;
    data.A = &A;
    data.x_vec = &x_vec;
    data.y_vec = &y_vec;

    return do_svds_nonsym(A_mat_r, m_scalar_r, n_scalar_r,
                          k_scalar_r, nu_scalar_r, nv_scalar_r,
                          params_list_r, sparse_mat_v_prod, sparse_mat_t_v_prod,
                          &data);

END_RCPP
}

