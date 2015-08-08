#ifndef MATTYPES_H
#define MATTYPES_H

#include "MatProd_matrix.h"

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
