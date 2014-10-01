#include <RcppEigen.h>
#include "EigsSym.h"
#include "EigsGen.h"
#include "MatTypes.h"

using Rcpp::as;
using std::string;

RcppExport SEXP eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP lower_logical_r,
                         SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    
    int n = as<int>(n_scalar_r);
    int nev = as<int>(k_scalar_r);
    int ncv = as<int>(params_rcpp["ncv"]);
    string which = as<string>(params_rcpp["which"]);
    int workmode = as<int>(params_rcpp["workmode"]);
    double sigma = as<double>(params_rcpp["sigma"]);
    char bmat = 'I';
    double tol = as<double>(params_rcpp["tol"]);
    int maxitr = as<int>(params_rcpp["maxitr"]);
    char uplo = as<bool>(lower_logical_r) ? 'L' : 'U';
    bool needSolve = (workmode != 1);
    bool retvec = as<bool>(params_rcpp["retvec"]);

    MatOp *op;
    switch(as<int>(mattype_scalar_r))
    {
        case (int) MATRIX:
            op = new MatOp_symmatrix(A_mat_r, n, uplo, sigma, needSolve);
            break;
        case (int) DSYMATRIX:
            op = new MatOp_dsyMatrix(A_mat_r, n, uplo, sigma, needSolve);
            break;
        default:
            Rcpp::stop("unsupported matrix type in eigs_sym()");
    }

    EigsSym eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.update();
    SEXP res = eig.extract(retvec);

    delete op;

    return res;

    END_RCPP
}


RcppExport SEXP eigs_gen(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    
    int n = as<int>(n_scalar_r);
    int nev = as<int>(k_scalar_r);
    int ncv = as<int>(params_rcpp["ncv"]);
    string which = as<string>(params_rcpp["which"]);
    int workmode = as<int>(params_rcpp["workmode"]);
    double sigmar = as<double>(params_rcpp["sigmar"]);
    double sigmai = as<double>(params_rcpp["sigmai"]);
    char bmat = 'I';
    double tol = as<double>(params_rcpp["tol"]);
    int maxitr = as<int>(params_rcpp["maxitr"]);
    bool needSolve = (workmode != 1);
    bool retvec = as<bool>(params_rcpp["retvec"]);

    MatOp *op;
    switch(as<int>(mattype_scalar_r))
    {
        case (int) MATRIX:
            op = new MatOp_matrix(A_mat_r, n, n, sigmar, sigmai, needSolve);
            break;
        case (int) DGEMATRIX:
            op = new MatOp_dgeMatrix(A_mat_r, n, n, sigmar, sigmai, needSolve);
            break;
        case (int) DGCMATRIX:
            op = new MatOp_dgCMatrix(A_mat_r, n, n, sigmar, sigmai, needSolve);
            break;
        case (int) DGRMATRIX:
            op = new MatOp_dgRMatrix(A_mat_r, n, n, sigmar, sigmai, needSolve);
            break;
        default:
            Rcpp::stop("unsupported matrix type in eigs()");
    }

    EigsGen eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.update();
    SEXP res = eig.extract(retvec);

    delete op;

    return res;

    END_RCPP
}
