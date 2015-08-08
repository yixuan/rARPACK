#ifndef MATTYPES_H
#define MATTYPES_H

#include "MatProd_matrix.h"
#include "MatProd_symmatrix.h"
#include "MatProd_dgeMatrix.h"
#include "MatProd_dsyMatrix.h"
#include "MatProd_sparseMatrix.h"

#include "RealShift_matrix.h"

enum MAT_TYPE {
    MATRIX = 0,
    SYMMATRIX,
    DGEMATRIX,
    DSYMATRIX,
    DGCMATRIX,
    DGRMATRIX,
    FUNCTION
};

#endif // MATTYPES_H
