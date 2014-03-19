#include <RcppEigen.h>
#include "EigsSymDense.h"

using std::string;
using Rcpp::as;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::PartialPivLU;

typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;

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

    EigsSymDense eig(n, nev, ncv, which, workmode,
                     sigma, bmat, tol, maxitr);
    eig.BindMatrix(A_mat_r, uplo);
    eig.Update();

    return eig.Extract();

    END_RCPP
}

