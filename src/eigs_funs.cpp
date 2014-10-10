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

    MatOp *op = newMatOp(A_mat_r, as<int>(mattype_scalar_r),
                         n, n, sigma, 0.0, true, needSolve,
                         uplo);

    EigsSym eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.compute(retvec);
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

    MatOp *op = newMatOp(A_mat_r, as<int>(mattype_scalar_r),
                         n, n, sigmar, sigmai, true, needSolve);

    EigsGen eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.compute(retvec);
    SEXP res = eig.extract();

    delete op;

    return res;

    END_RCPP
}

RcppExport SEXP eigs_fun(SEXP FUN_function_r, SEXP args_list_r,
                         SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);
    
    int n = as<int>(n_scalar_r);
    int nev = as<int>(k_scalar_r);
    int ncv = as<int>(params_rcpp["ncv"]);
    string which = as<string>(params_rcpp["which"]);
    int workmode = as<int>(params_rcpp["workmode"]);
    // To be safe, we force workmode to be 1
    workmode = 1;
    char bmat = 'I';
    double tol = as<double>(params_rcpp["tol"]);
    int maxitr = as<int>(params_rcpp["maxitr"]);
    bool retvec = as<bool>(params_rcpp["retvec"]);

    MatOp *op = newMatOp(FUN_function_r, as<int>(mattype_scalar_r),
                         n, n, 0.0, 0.0, false, false, '\0', args_list_r);

    EigsGen eig(n, nev, ncv, op, which, workmode,
                bmat, tol, maxitr);
    eig.compute(retvec);
    SEXP res = eig.extract();

    delete op;

    return res;

    END_RCPP
}
