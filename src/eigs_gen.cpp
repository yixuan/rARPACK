#include <RcppEigen.h>
#include <GenEigsSolver.h>
#include "MatOp/MatTypes.h"

using Rcpp::as;

enum SOLVER_TYPE {
    REGULAR = 0,
    REAL_SHIFT,
    COMPLEX_SHIFT
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
GenEigsSolver<double, RULE, OpType> eigs(&op, nev, ncv);                       \
EIG_COMMON_CODE



#define EIG_CODE_REAL_SHIFT(RULE)                                              \
GenEigsRealShiftSolver<double, RULE, OpType> eigs(&op, nev, ncv, sigmar);      \
EIG_COMMON_CODE



#define EIG_CODE_COMPLEX_SHIFT(RULE)                                           \
GenEigsComplexShiftSolver<double, RULE, OpType> eigs(&op, nev, ncv, sigmar, sigmai);      \
EIG_COMMON_CODE



#define EIG_CODE_GENERATOR(SOLVER)                                             \
switch(rule)                                                                   \
{                                                                              \
    case WHICH_LM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LM) }                                      \
        break;                                                                 \
    case WHICH_LR :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LR) }                                      \
        break;                                                                 \
    case WHICH_LI :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LI) }                                      \
        break;                                                                 \
    case WHICH_SM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SM) }                                      \
        break;                                                                 \
    case WHICH_SR :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SR) }                                      \
        break;                                                                 \
    case WHICH_SI :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SI) }                                      \
        break;                                                                 \
    default:                                                                   \
        Rcpp::stop("unsupported selection rule");                              \
}

/************************ Macros to generate code ************************/



/************************ Regular mode ************************/
template <typename OpType>
Rcpp::RObject run_eigs_gen(OpType &op, const int rule, const double *init_resid,
                           int nev, int ncv, int maxitr, double tol, bool retvec)
{
    int nconv;
    Rcpp::RObject evals, evecs;

    EIG_CODE_GENERATOR(REGULAR)

    return R_NilValue;  // should not reach here
}


RcppExport SEXP eigs_gen(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                         SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n        = as<int>(n_scalar_r);
    int nev      = as<int>(k_scalar_r);
    int ncv      = as<int>(params_rcpp["ncv"]);
    int rule     = as<int>(params_rcpp["which"]);
    double tol   = as<double>(params_rcpp["tol"]);
    int maxitr   = as<int>(params_rcpp["maxitr"]);
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
            res = run_eigs_gen(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
        }
        break;
        case DGEMATRIX:
        {
            MatProd_dgeMatrix op(A_mat_r, n, n);
            res = run_eigs_gen(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
        }
        break;
        case DGCMATRIX:
        {
            MatProd_dgCMatrix op(A_mat_r, n, n);
            res = run_eigs_gen(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
        }
        break;
        case DGRMATRIX:
        {
            MatProd_dgRMatrix op(A_mat_r, n, n);
            res = run_eigs_gen(op, rule, init_resid, nev, ncv, maxitr, tol, retvec);
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
/************************ Regular mode ************************/
