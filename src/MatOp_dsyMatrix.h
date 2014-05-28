#ifndef MATOPDSYMATRIX_H
#define MATOPDSYMATRIX_H

#include "MatOp_symmatrix.h"

inline char get_uplo(SEXP dsyMatrix)
{
    SEXP str = GET_SLOT(dsyMatrix, Rf_install("uplo"));
    return CHAR(STRING_ELT(str, 0))[0];
}

// Operations on "dsyMatrix" class, defined in Matrix package
class MatOp_dsyMatrix : public MatOp_symmatrix
{
public:
    // Constructor
    MatOp_dsyMatrix(SEXP mat_, double sigma_ = 0,
                    bool needSolve_ = false):
        MatOp_symmatrix(GET_SLOT(mat_, Rf_install("x")),
                        get_uplo(mat_),
                        sigma_,
                        needSolve_)
    {}
    virtual ~MatOp_dsyMatrix() {}
};


#endif // MATOPDSYMATRIX_H
