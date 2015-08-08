#include <RcppEigen.h>
#include <SymEigsSolver.h>
#include "MatOp/MatTypes.h"

using Rcpp::as;
using std::string;

enum SOLVER_TYPE {
    REGULAR = 0,
    REAL_SHIFT
};

/************************ Macros to generate code ************************/

#define EIG_COMMON_CODE                                                        \
eigs.init(init_resid);                                                         \
nconv = eigs.compute(maxitr, tol);                                             \
if(nconv < nev)                                                                \
    Rcpp::warning("only %d eigenvalues converged, less than k", nconv);        \
evals = Rcpp::wrap(eigs.eigenvalues());                                        \
if(retvec)                                                                     \
    evecs = Rcpp::wrap(eigs.eigenvectors());                                   \
else                                                                           \
    evecs = R_NilValue;                                                        \
                                                                               \
return Rcpp::List::create(                                                     \
    Rcpp::Named("values")  = evals,                                            \
    Rcpp::Named("vectors") = evecs,                                            \
    Rcpp::Named("nconv")   = nconv,                                            \
    Rcpp::Named("niter")   = eigs.num_iterations(),                            \
    Rcpp::Named("nops")    = eigs.num_operations()                             \
);



#define EIG_CODE_REGULAR(RULE)                                                 \
SymEigsSolver<double, RULE, OpType> eigs(&op, nev, ncv);                       \
EIG_COMMON_CODE



#define EIG_CODE_REAL_SHIFT(RULE)                                              \
SymEigsShiftSolver<double, RULE, OpType> eigs(&op, nev, ncv, sigma);           \
EIG_COMMON_CODE



#define EIG_CODE_GENERATOR(SOLVER)                                             \
switch(rule)                                                                   \
{                                                                              \
    case WHICH_LM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LM) }                                      \
        break;                                                                 \
    case WHICH_LA :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LA) }                                      \
        break;                                                                 \
    case WHICH_SM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SM) }                                      \
        break;                                                                 \
    case WHICH_SA :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SA) }                                      \
        break;                                                                 \
    case WHICH_BE :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_BE) }                                      \
        break;                                                                 \
    default:                                                                   \
        Rcpp::stop("unsupported selection rule");                              \
}

/************************ Macros to generate code ************************/

template <typename OpType>
Rcpp::RObject run_eigs_sym(OpType &op, const int rule, const double *init_resid,
                           int nev, int ncv, int maxitr, double tol, bool retvec)
{
    int nconv;
    Rcpp::RObject evals, evecs;

    EIG_CODE_GENERATOR(REGULAR)

    return R_NilValue;  // should not reach here
}


RcppExport SEXP eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP lower_logical_r,
                         SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n        = as<int>(n_scalar_r);
    int nev      = as<int>(k_scalar_r);
    int ncv      = as<int>(params_rcpp["ncv"]);
    int rule     = as<int>(params_rcpp["which"]);
    double tol   = as<double>(params_rcpp["tol"]);
    int maxitr   = as<int>(params_rcpp["maxitr"]);
    char uplo    = as<bool>(lower_logical_r) ? 'L' : 'U';
    bool retvec  = as<bool>(params_rcpp["retvec"]);
    int mattype  = as<int>(mattype_scalar_r);

    // Prepare initial residuals
    double *init_resid;
    #include "rands.h"
    if(n <= rands_len)
    {
        init_resid = rands;
    } else {
        init_resid = new double[n];
        double *coef_pntr = init_resid;
        for(int i = 0; i < n / rands_len; i++, coef_pntr += rands_len)
        {
            std::copy(rands, rands + rands_len, coef_pntr);
        }
        std::copy(rands, rands + n % rands_len, coef_pntr);
    }

    Rcpp::RObject res;

    switch(mattype)
    {
        case MATRIX:
        {
            MatProd_matrix op(A_mat_r, n, n);
            res = run_eigs_sym(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
        }
        break;
        case SYMMATRIX:
        {
            MatProd_symmatrix op(A_mat_r, n, uplo);
            res = run_eigs_sym(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
        }
        break;
        case DGEMATRIX:
        {
            MatProd_dgeMatrix op(A_mat_r, n, n);
            res = run_eigs_sym(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
        }
        break;
        default:
            Rcpp::stop("unsupported matrix type");
    }

    if(n > rands_len)
        delete [] init_resid;

    return res;  // should not reach here

    END_RCPP
}
