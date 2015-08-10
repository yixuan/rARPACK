#include <RcppEigen.h>
#include <GenEigsSolver.h>
#include <MatOp/MatTypes.h>

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
evals = Rcpp::wrap(eigs.eigenvalues());                                        \
if(retvec)                                                                     \
    evecs = Rcpp::wrap(eigs.eigenvectors());                                   \
else                                                                           \
    evecs = R_NilValue;                                                        \
niter = eigs.num_iterations();                                                 \
nops = eigs.num_operations();



#define EIG_CODE_REGULAR(RULE, OPTYPE)                                         \
GenEigsSolver<double, RULE, OPTYPE> eigs(op, nev, ncv);                        \
EIG_COMMON_CODE



#define EIG_CODE_REAL_SHIFT(RULE, OPTYPE)                                      \
GenEigsRealShiftSolver<double, RULE, OPTYPE> eigs(op, nev, ncv, sigmar);       \
EIG_COMMON_CODE



#define EIG_CODE_COMPLEX_SHIFT(RULE, OPTYPE)                                   \
GenEigsComplexShiftSolver<double, RULE, OPTYPE> eigs(op, nev, ncv, sigmar, sigmai);      \
EIG_COMMON_CODE



#define EIG_CODE_GENERATOR(SOLVER, OPTYPE)                                     \
switch(rule)                                                                   \
{                                                                              \
    case WHICH_LM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LM, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_LR :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LR, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_LI :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_LI, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SM :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SM, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SR :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SR, OPTYPE) }                              \
        break;                                                                 \
    case WHICH_SI :                                                            \
        { EIG_CODE_ ## SOLVER(WHICH_SI, OPTYPE) }                              \
        break;                                                                 \
    default:                                                                   \
        Rcpp::stop("unsupported selection rule");                              \
}

/************************ Macros to generate code ************************/



/************************ Regular mode ************************/
Rcpp::RObject run_eigs_gen(MatProd* op, int n, int nev, int ncv, int rule,
                           int maxitr, double tol, bool retvec)
{
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

    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(REGULAR, MatProd)

    if(nconv < nev)
        Rcpp::warning("only %d eigenvalues converged, less than k", nconv);

    if(n > rands_len)
        delete [] init_resid;

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
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

    MatProd *op;
    Rcpp::RObject res;

    switch(mattype)
    {
        case MATRIX:
            op = new MatProd_matrix(A_mat_r, n, n);
            break;
        case DGEMATRIX:
            op = new MatProd_dgeMatrix(A_mat_r, n, n);
            break;
        case DGCMATRIX:
            op = new MatProd_dgCMatrix(A_mat_r, n, n);
            break;
        case DGRMATRIX:
            op = new MatProd_dgRMatrix(A_mat_r, n, n);
            break;
        default:
            Rcpp::stop("unsupported matrix type");
            // Eliminate compiler warning, but should not reach here
            op = new MatProd_matrix(A_mat_r, n, n);
    }

    res = run_eigs_gen(op, n, nev, ncv, rule, maxitr, tol, retvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Regular mode ************************/



/************************ Real shift mode ************************/
Rcpp::RObject run_eigs_real_shift_gen(RealShift* op, int n, int nev, int ncv, int rule,
                                      double sigmar, int maxitr, double tol, bool retvec)
{
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

    Rcpp::RObject evals, evecs;
    int nconv = 0, niter = 0, nops = 0;

    EIG_CODE_GENERATOR(REAL_SHIFT, RealShift)

    if(nconv < nev)
        Rcpp::warning("only %d eigenvalues converged, less than k", nconv);

    if(n > rands_len)
        delete [] init_resid;

    return Rcpp::List::create(
        Rcpp::Named("values")  = evals,
        Rcpp::Named("vectors") = evecs,
        Rcpp::Named("nconv")   = nconv,
        Rcpp::Named("niter")   = niter,
        Rcpp::Named("nops")    = nops
    );
}

RcppExport SEXP eigs_real_shift_gen(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                                    SEXP params_list_r, SEXP mattype_scalar_r)
{
    BEGIN_RCPP

    Rcpp::List params_rcpp(params_list_r);

    int n         = as<int>(n_scalar_r);
    int nev       = as<int>(k_scalar_r);
    int ncv       = as<int>(params_rcpp["ncv"]);
    int rule      = as<int>(params_rcpp["which"]);
    double tol    = as<double>(params_rcpp["tol"]);
    int maxitr    = as<int>(params_rcpp["maxitr"]);
    bool retvec   = as<bool>(params_rcpp["retvec"]);
    int mattype   = as<int>(mattype_scalar_r);
    double sigmar = as<double>(params_rcpp["sigmar"]);

    RealShift *op;
    Rcpp::RObject res;

    switch(mattype)
    {
        case MATRIX:
            op = new RealShift_matrix(A_mat_r, n);
            break;
        case DGEMATRIX:
            op = new RealShift_dgeMatrix(A_mat_r, n);
            break;
        case DGCMATRIX:
            op = new RealShift_dgCMatrix(A_mat_r, n);
            break;
        case DGRMATRIX:
            op = new RealShift_dgRMatrix(A_mat_r, n);
            break;
        default:
            Rcpp::stop("unsupported matrix type");
            // Eliminate compiler warning, but should not reach here
            op = new RealShift_matrix(A_mat_r, n);
    }

    res = run_eigs_real_shift_gen(op, n, nev, ncv, rule, sigmar, maxitr, tol, retvec);

    delete op;

    return res;

    END_RCPP
}
/************************ Real shift mode ************************/
