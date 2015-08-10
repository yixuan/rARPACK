#ifndef MATPROD_DSYMATRIX_H
#define MATPROD_DSYMATRIX_H

#include <Rcpp.h>
#include <Rdefines.h>  // for R macros
#include "MatProd_symmatrix.h"

class MatProd_dsyMatrix : public MatProd_symmatrix
{
public:
    MatProd_dsyMatrix(SEXP mat_, const int nrow_, const char uplo_ = 'L') :
        MatProd_symmatrix(GET_SLOT(mat_, Rf_install("x")), nrow_, uplo_)
    {}
};


#endif // MATPROD_DSYMATRIX_H
