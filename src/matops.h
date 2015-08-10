#include <rarpack/MatOps.h>
#include <Rinternals.h>

enum MAT_TYPE {
    MATRIX = 0,
    SYMMATRIX,
    DGEMATRIX,
    DSYMATRIX,
    DGCMATRIX,
    DGRMATRIX,
    FUNCTION
};

MatProd* get_mat_prod_op(SEXP mat, int n, SEXP extra_arg, int mat_type);

RealShift* eigs_sym_get_real_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type);

RealShift* eigs_gen_get_real_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type);

ComplexShift* get_complex_shift_op(SEXP mat, int n, SEXP extra_arg, int mat_type);
