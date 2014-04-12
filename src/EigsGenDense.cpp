#include "EigsGenDense.h"

const char EigsGenDense::BLAS_trans = 'N';
const double EigsGenDense::BLAS_alpha = 1.0;
const int EigsGenDense::BLAS_one = 1;
const double EigsGenDense::BLAS_zero = 0.0;

EigsGenDense::EigsGenDense(int n_, int nev_, int ncv_,
                           const string & which_, int workmode_,
                           double sigmar_, double sigmai_,
                           char bmat_, double tol_, int maxitr_) :
    EigsGen(n_, nev_, ncv_, which_, workmode_,
            sigmar_, sigmai_, bmat_, tol_, maxitr_),
    cx_vec(n),
    x_vec(NULL, n),
    y_vec(NULL, n)
{
    A_pntr = NULL;
}

EigsGenDense::~EigsGenDense()
{
    delete [] workv;
    delete [] workl;
    delete [] workd;
    delete [] resid;
}

void EigsGenDense::BindMatrix(SEXP mat_)
{
    A_pntr = REAL(mat_);

    // If we don't need shift-and-invert mode
    if (workmode == 1)
    {
        MatrixLinked(true);
        return;
    }

    // Map mat_ to Eigen matrix
    MapMat A(A_pntr, n, n);

    // If sigma is real
    if(fabs(sigmai) < 1e-17)
    {
        // Subtract the diagonal elements by sigma, i.e., A - sigma * I
        for(int i = 0; i < n; i++)
        {
            A(i, i) -= sigmar;
        }
        // LU decomposition
        solver.compute(A);
        // Change A back
        for(int i = 0; i < n; i++)
        {
            A(i, i) += sigmar;
        }
    } else {
        MatrixXcd cA = A.cast< std::complex<double> >();
        // Subtract the diagonal elements by sigma, i.e., A - sigma * I
        for(int i = 0; i < n; i++)
        {
            cA(i, i) -= std::complex<double>(sigmar, sigmai);
        }
        // LU decomposition
        csolver.compute(cA);
        cx_vec.setZero();
    }
    
    MatrixLinked(true);
}

void EigsGenDense::MultVector(double *x_in, double *y_out)
{
    F77_CALL(dgemv)(&BLAS_trans, &n, &n,
            &BLAS_alpha, A_pntr, &n,
            x_in, &BLAS_one, &BLAS_zero,
            y_out, &BLAS_one);
}

void EigsGenDense::MultVectorShift(double *x_in, double *y_out)
{
    if(fabs(sigmai) < 1e-17)
    {
        new (&x_vec) MapVec(x_in, n);
        new (&y_vec) MapVec(y_out, n);
        y_vec = solver.solve(x_vec);
    } else {
        cx_vec.real() = MapVec(x_in, n);
        new (&y_vec) MapVec(y_out, n);
        y_vec = csolver.solve(cx_vec).real();
    }
}

