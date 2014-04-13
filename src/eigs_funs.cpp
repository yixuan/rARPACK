#include <RcppEigen.h>
#include "EigsSym.h"
#include "EigsGen.h"
#include "MatOpMatrix.h"
#include "MatOpSymMatrix.h"
#include "MatOpDGCMatrix.h"

using Rcpp::as;

RcppExport SEXP den_real_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                             SEXP params_list_r,
                             SEXP lower_logical_r)
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

    MatOpSymMatrix op(A_mat_r, uplo, sigma, needSolve);
    EigsSym eig(n, nev, ncv, &op, which, workmode,
                bmat, tol, maxitr);
    eig.update();

    return eig.extract(retvec);

    END_RCPP
}

RcppExport SEXP den_real_gen(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                             SEXP params_list_r)
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

    MatOpMatrix op(A_mat_r, sigmar, sigmai, needSolve);
    EigsGen eig(n, nev, ncv, &op, which, workmode,
                bmat, tol, maxitr);
    eig.update();

    return eig.extract(retvec);

    END_RCPP
}

RcppExport SEXP sparse_real_gen(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                SEXP params_list_r)
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

    MatOpDGCMatrix op(A_mat_r, sigmar, sigmai, needSolve);
    EigsGen eig(n, nev, ncv, &op, which, workmode,
                bmat, tol, maxitr);
    eig.update();

    return eig.extract(retvec);

    END_RCPP
}
