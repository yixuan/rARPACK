#include <RcppEigen.h>
#include <Rdefines.h>
#include "EigsSym.h"
#include "EigsGen.h"
#include "MatTypes.h"

using Rcpp::as;

RcppExport SEXP eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP lower_logical_r,
                         SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    
    int n = INTEGER(n_scalar_r)[0];
    int nev = INTEGER(k_scalar_r)[0];
    int ncv = as<int>(params_rcpp["ncv"]);
    Rcpp::CharacterVector which_rcpp = params_rcpp["which"];
    string which = string(which_rcpp[0]);
    int workmode = as<int>(params_rcpp["workmode"]);
    double sigma = as<double>(params_rcpp["sigma"]);
    char bmat = 'I';
    double tol = as<double>(params_rcpp["tol"]);
    int maxitr = as<int>(params_rcpp["maxitr"]);
    char uplo = LOGICAL(lower_logical_r)[0] ? 'L' : 'U';
    bool needSolve = (workmode != 1);
    bool retvec = as<bool>(params_rcpp["retvec"]);

    MatOp *op;
    switch(INTEGER(mattype_scalar_r)[0])
    {
        case (int) MATRIX:
            op = new MatOp_symmatrix(A_mat_r, uplo, sigma, needSolve);
            break;
        case (int) DSYMATRIX:
            op = new MatOp_dsyMatrix(A_mat_r, sigma, needSolve);
            break;
        default:
            Rcpp::stop("unsupported matrix type in eigs_sym()");
    }

    EigsSym eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.update();
    SEXP res = eig.extract();

    delete op;

    return res;

    END_RCPP
}


RcppExport SEXP eigs_gen(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    
    int n = INTEGER(n_scalar_r)[0];
    int nev = INTEGER(k_scalar_r)[0];
    int ncv = as<int>(params_rcpp["ncv"]);
    Rcpp::CharacterVector which_rcpp = params_rcpp["which"];
    string which = string(which_rcpp[0]);
    int workmode = as<int>(params_rcpp["workmode"]);
    double sigmar = as<double>(params_rcpp["sigmar"]);
    double sigmai = as<double>(params_rcpp["sigmai"]);
    char bmat = 'I';
    double tol = as<double>(params_rcpp["tol"]);
    int maxitr = as<int>(params_rcpp["maxitr"]);
    bool needSolve = (workmode != 1);
    bool retvec = as<bool>(params_rcpp["retvec"]);

    MatOp *op;
    switch(INTEGER(mattype_scalar_r)[0])
    {
        case (int) MATRIX:
            op = new MatOp_matrix(A_mat_r, sigmar, sigmai, needSolve);
            break;
        case (int) DGEMATRIX:
            op = new MatOp_dgeMatrix(A_mat_r, sigmar, sigmai, needSolve);
            break;
        case (int) DGCMATRIX:
            op = new MatOp_dgCMatrix(A_mat_r, sigmar, sigmai, needSolve);
            break;
        default:
            Rcpp::stop("unsupported matrix type in eigs()");
    }

    EigsGen eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.update();
    SEXP res = eig.extract();

    delete op;

    return res;

    END_RCPP
}
