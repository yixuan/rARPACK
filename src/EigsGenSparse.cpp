#include "EigsGenSparse.h"

EigsGenSparse::EigsGenSparse(int n_, int nev_, int ncv_,
                             const string & which_, int workmode_,
                             double sigmar_, double sigmai_,
                             char bmat_, double tol_, int maxitr_) :
    EigsGen(n_, nev_, ncv_, which_, workmode_,
            sigmar_, sigmai_, bmat_, tol_, maxitr_),
    A(n, n, 0, NULL, NULL, NULL),
    x_vec(NULL, n),
    y_vec(NULL, n)
{
}

EigsGenSparse::~EigsGenSparse()
{
    delete [] workv;
    delete [] workl;
    delete [] workd;
    delete [] resid;
}

void EigsGenSparse::BindMatrix(SEXP mat_)
{
    A = Rcpp::as<MapSpMat>(mat_);

    // If we don't need shift-and-invert mode
    if (workmode == 1)
    {
        MatrixLinked(true);
        return;
    }

    // If sigma is real
    if(fabs(sigmai) < 1e-17)
    {
        // Create a sparse idendity matrix
        SpMat I(n, n);
        I.setIdentity();

        // Sparse LU decomposition
        solver.compute(A - sigmar * I);
        if(solver.info() != Eigen::Success) {
            Rcpp::stop("Eigen: LU decomposition of A - sigma * I failed");
        }

    } else {
        SpCMat cA = A.cast< std::complex<double> >();
        
        // Create a sparse identity matrix (1 + 0i on diagonal)
        SpCMat I(n, n);
        I.setIdentity();
        
        // Sparse LU decomposition
        csolver.compute(cA - std::complex<double>(sigmar, sigmai) * I);
        if(csolver.info() != Eigen::Success) {
            Rcpp::stop("Eigen: LU decomposition of A - sigma * I failed");
        }

        cx_vec.resize(n);
        cx_vec.setZero();
    }
    
    MatrixLinked(true);
}

void EigsGenSparse::MultVector(double *x_in, double *y_out)
{
    new (&x_vec) MapVec(x_in, n);
    new (&y_vec) MapVec(y_out, n);
    y_vec = A * x_vec;
}

void EigsGenSparse::MultVectorShift(double *x_in, double *y_out)
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

