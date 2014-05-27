#include <RcppEigen.h>
#include <Rdefines.h>
#include "ARPACK.h"
#include "EigsSym.h"
#include "MatOp_matrix.h"
#include "MatOp_symmatrix.h"
#include "MatOp_dsyMatrix.h"
#include "MatOp_dgCMatrix.h"

SEXP do_svds_sym(MatOp *op, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r);

SEXP do_svds_gen(MatOp *op, SEXP m_scalar_r, SEXP n_scalar_r,
                 SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                 SEXP params_list_r);

enum SVDSSYMMAT {
    DSYMATRIX = 0
};

enum SVDSGENMAT {
    MATRIX = 0,
    DGCMATRIX
};


// Main function to calculate truncated SVD for
// dense, real, symmetric matrices
RcppExport SEXP svds_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
BEGIN_RCPP

    MatOp *op;
    switch(INTEGER(mattype_scalar_r)[0])
    {
        case (int) DSYMATRIX:
            op = new MatOp_dsyMatrix(A_mat_r, 0, false);
            break;
        default:
            Rcpp::stop("unsupported matrix type in svds()");
    }

    SEXP res = do_svds_sym(op, n_scalar_r, k_scalar_r,
                           nu_scalar_r, nv_scalar_r,
                           params_list_r);

    delete op;

    return res;

END_RCPP
}

// Main function to calculate truncated SVD for
// real general matrices
RcppExport SEXP svds_gen(SEXP A_mat_r, SEXP m_scalar_r, SEXP n_scalar_r,
                         SEXP k_scalar_r, SEXP nu_scalar_r, SEXP nv_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
BEGIN_RCPP

    MatOp *op;
    switch(INTEGER(mattype_scalar_r)[0])
    {
        case (int) MATRIX:
            op = new MatOp_matrix(A_mat_r, 0, 0, false);
            break;
        case (int) DGCMATRIX:
            op = new MatOp_dgCMatrix(A_mat_r, 0, 0, false);
            break;
        default:
            Rcpp::stop("unsupported matrix type in svds()");
    }

    SEXP res = do_svds_gen(op, m_scalar_r, n_scalar_r, k_scalar_r,
                           nu_scalar_r, nv_scalar_r,
                           params_list_r);

    delete op;

    return res;

END_RCPP
}
