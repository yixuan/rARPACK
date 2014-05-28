#ifndef MATOPDGEMATRIX_H
#define MATOPDGEMATRIX_H

#include "MatOp_matrix.h"

// Operations on "dgeMatrix" class, defined in Matrix package
class MatOp_dgeMatrix : public MatOp_matrix
{
public:
    // Constructor
    MatOp_dgeMatrix(SEXP mat_, double sigmar_ = 0, double sigmai_ = 0,
                    bool needSolve_ = false) :
        MatOp_matrix(REAL(GET_SLOT(mat_, Rf_install("x"))),
                     INTEGER(GET_SLOT(mat_, Rf_install("Dim")))[0],
                     INTEGER(GET_SLOT(mat_, Rf_install("Dim")))[1],
                     sigmar_, sigmai_, needSolve_)
    {}
    virtual ~MatOp_dgeMatrix() {}
};


#endif // MATOPDGEMATRIX_H
