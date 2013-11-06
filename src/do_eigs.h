#ifndef DO_EIGS_H
#define DO_EIGS_H

#include <Rcpp.h>
#include <R_ext/BLAS.h>

#include "arpack_f.h"
#include "arpack_cpp.h"

// Prototype of function to calculate
// matrix(n, n) * vector(n) product
// A: an n*n matrix, may be sparse
// x_in: an array of length n
// y_out: result of A * x_in
// data: additional data passed to function
typedef void (*Mvfun)(SEXP mat, double *x_in, double *y_out,
                      int n, void *data);

// Common function to calculate eigen values/vectors
// mat_v_prod and data should be implemented according to
// the problem (whether matrix is dense of sparse, for example)
SEXP do_eigs_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                    SEXP params_list_r,
                    Mvfun mat_v_prod, void *data);

SEXP do_eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                 SEXP params_list_r,
                 Mvfun mat_v_prod, void *data);


#endif // DO_EIGS_H
