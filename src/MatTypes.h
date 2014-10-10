#ifndef MATTYPES_H
#define MATTYPES_H

#include "MatOp.h"
#include "MatOp_matrix.h"
#include "MatOp_symmatrix.h"
#include "MatOp_dgeMatrix.h"
#include "MatOp_dsyMatrix.h"
#include "MatOp_sparseMatrix.h"
#include "MatOp_function.h"

enum MATTYPE {
    MATRIX = 0,
    SYMMATRIX,
    DGEMATRIX,
    DSYMATRIX,
    DGCMATRIX,
    DGRMATRIX,
    FUNCTION
};

inline MatOp* newMatOp(SEXP mat, int mattype, int m, int n,
                       double sigmar = 0.0, double sigmai = 0.0,
                       bool needTprod = true, bool needSolve = false,
                       char uplo = 'L', SEXP args = R_NilValue)
{
    MatOp *op = NULL;
    switch(mattype)
    {
        case (int) MATRIX:
            op = new MatOp_matrix(mat, m, n, sigmar, sigmai, needSolve);
            break;
        case (int) SYMMATRIX:
            op = new MatOp_symmatrix(mat, n, uplo, sigmar, needSolve);
            break;
        case (int) DGEMATRIX:
            op = new MatOp_dgeMatrix(mat, m, n, sigmar, sigmai, needSolve);
            break;
        case (int) DSYMATRIX:
            op = new MatOp_dsyMatrix(mat, n, uplo, sigmar, needSolve);
            break;
        case (int) DGCMATRIX:
            op = new MatOp_dgCMatrix(mat, m, n, sigmar, sigmai, needSolve);
            break;
        case (int) DGRMATRIX:
            op = new MatOp_dgRMatrix(mat, m, n, sigmar, sigmai, needSolve);
            break;
        case (int) FUNCTION:
            op = new MatOp_function(mat, n, args);
            break;
        default:
            Rcpp::stop("unsupported matrix type");
    }
    return op;
}

#endif // MATTYPES_H
