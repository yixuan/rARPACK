#ifndef MATTYPES_H
#define MATTYPES_H

#include "MatOp/MatProd.h"
#include "MatOp/MatProd_matrix.h"
#include "MatOp/MatProd_symmatrix.h"
#include "MatOp/MatProd_dgeMatrix.h"
#include "MatOp/MatProd_dsyMatrix.h"
#include "MatOp/MatProd_sparseMatrix.h"

#include "MatOp/RealShift.h"
#include "MatOp/RealShift_matrix.h"
#include "MatOp/RealShift_symmatrix.h"
#include "MatOp/RealShift_dgeMatrix.h"
#include "MatOp/RealShift_dsyMatrix.h"
#include "MatOp/RealShift_sparseMatrix.h"

#include "MatOp/ComplexShift.h"
#include "MatOp/ComplexShift_matrix.h"
#include "MatOp/ComplexShift_dgeMatrix.h"
#include "MatOp/ComplexShift_sparseMatrix.h"

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
