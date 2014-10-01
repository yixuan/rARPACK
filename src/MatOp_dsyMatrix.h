#ifndef MATOP_DSYMATRIX_H
#define MATOP_DSYMATRIX_H

#include "MatOp_symmatrix.h"
#include <Rdefines.h>  // for R macros

// Operations on "dsyMatrix" class, defined in Matrix package
class MatOp_dsyMatrix : public MatOp_symmatrix
{
public:
    // Constructor
    MatOp_dsyMatrix(SEXP mat_, int n_, char uplo_ = 'L',
                    double sigma_ = 0,
                    bool needSolve_ = false):
        MatOp_symmatrix(GET_SLOT(mat_, Rf_install("x")),
                        n_, uplo_, sigma_, needSolve_)
    {
        if((uplo != 'L') && (uplo != 'U'))
            uplo = get_uplo(mat_);
    }
    virtual ~MatOp_dsyMatrix() {}
    
    static char get_uplo(SEXP dsyMatrix)
    {
        SEXP str = GET_SLOT(dsyMatrix, Rf_install("uplo"));
        return CHAR(STRING_ELT(str, 0))[0];
    }
};


#endif // MATOP_DSYMATRIX_H
