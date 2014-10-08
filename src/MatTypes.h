#ifndef MATTYPES_H
#define MATTYPES_H

#include "MatOp.h"
#include "MatOp_matrix.h"
#include "MatOp_symmatrix.h"
#include "MatOp_dgeMatrix.h"
#include "MatOp_dsyMatrix.h"
#include "MatOp_sparseMatrix.h"
#include "MatOp_function.h"

enum MATTYPE {
    MATRIX = 0,
    DGEMATRIX,
    DSYMATRIX,
    DGCMATRIX,
    DGRMATRIX,
    FUNCTION
};

#endif // MATTYPES_H
