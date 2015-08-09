#ifndef MATTYPES_H
#define MATTYPES_H

#include "MatProd.h"
#include "MatProd_matrix.h"
#include "MatProd_symmatrix.h"
#include "MatProd_dgeMatrix.h"
#include "MatProd_dsyMatrix.h"
#include "MatProd_sparseMatrix.h"

#include "RealShift_matrix.h"
#include "RealShift_symmatrix.h"
#include "RealShift_dgeMatrix.h"
#include "RealShift_dsyMatrix.h"
#include "RealShift_sparseMatrix.h"

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
