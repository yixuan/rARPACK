#ifndef MATOP_DGEMATRIX_H
#define MATOP_DGEMATRIX_H

#include "MatOp_matrix.h"
#include <Rdefines.h>  // for R macros

// Operations on "dgeMatrix" class, defined in Matrix package
class MatOp_dgeMatrix : public MatOp_matrix
{
public:
    // Constructor
    MatOp_dgeMatrix(SEXP mat_, int m_, int n_,
                    double sigmar_ = 0, double sigmai_ = 0,
                    bool needSolve_ = false) :
        MatOp_matrix(REAL(GET_SLOT(mat_, Rf_install("x"))),
                     m_, n_, sigmar_, sigmai_, needSolve_)
    {}
    virtual ~MatOp_dgeMatrix() {}
};


#endif // MATOP_DGEMATRIX_H
