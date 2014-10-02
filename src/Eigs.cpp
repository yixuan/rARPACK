#include "Eigs.h"
#include <utility>
#include <vector>
#include <algorithm>

using std::string;

Eigs::Eigs(int n_, int nev_, int ncv_, MatOp *op_,
           const string & which_, int workmode_,
           char bmat_, double tol_, int maxitr_) :
    n(n_), nev(nev_), ncv(ncv_), op(op_),
    which(which_), workmode(workmode_),
    bmat(bmat_), tol(tol_), maxitr(maxitr_),
    ido(0), info(0), ierr(0), retvec(false)
{
    for(int i = 0; i < 11; i++)
        iparam[i] = 0;

    iparam[1 - 1] = 1;
    iparam[3 - 1] = maxitr;
    iparam[7 - 1] = workmode;

    for(int i = 0; i < 14; i++)
        ipntr[i] = 0;

    resid = new double[n]();
    workd = new double[3 * n]();
    // lworkl and workl are not initialized here since their
    // sizes are determined according to the problem.
    // Classes derived from Eigs should do this.
}

unsigned long int Eigs::seed_next = 1;

void Eigs::srand(unsigned int seed)
{
    seed_next = seed;
}

unsigned int Eigs::rand()
{
    seed_next = seed_next * 1103515245 + 12345;
    return (unsigned int)(seed_next/65536) % 32768;
}

void Eigs::initResid()
{
    // Create initial residual vector
    // info = 1 means using the residual vector we provide
    info = 1;
    double *initcoef = new double[n];
    Rcpp::RNGScope scp;
    for(int i = 0; i < n; i++)
        initcoef[i] = R::unif_rand() - 0.5;

    // resid = A * initcoef
    op->prod(initcoef, resid);
    delete [] initcoef;
}

void Eigs::matOp()
{
    if (workmode == 3)  // Shift-and-invert mode
        op->shiftSolve(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
    else                // Regular mode
        op->prod(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
}

void Eigs::compute(bool rvec)
{
    initResid();

    while (ido != 99)
    {
        aupd();
        switch(ido)
        {
            case -1:
            case 1:
                matOp();
                break;
            default:
                break;
        }
    }
    
    // info > 0 means warning, < 0 means error
    if (info > 0)  warning(1, info);
    if (info < 0)  error(1, info);
    
    retvec = rvec;
    eupd();
    
    // ierr > 0 means warning, < 0 means error
    if (ierr > 0) warning(2, ierr);
    if (ierr < 0) error(2, ierr);
}

Eigs::~Eigs()
{
    delete [] workd;
    delete [] resid;
}



// Copy source[, i] to dest[, j]
void copy_column(const Rcpp::NumericMatrix &source, int i,
                 Rcpp::NumericMatrix &dest, int j)
{
    int n1 = source.nrow();
    int n2 = dest.nrow();
    if(n1 != n2)  return;
    
    std::copy(&source(0, i), &source(0, i) + n1, &dest(0, j));
}
// Copy src_real[, i] and src_img[, j] to dest[, k]
void copy_column(const Rcpp::NumericMatrix &src_real, int i,
                 const Rcpp::NumericMatrix &src_img, int j,
                 Rcpp::ComplexMatrix &dest_comp, int k)
{
    int n1r = src_real.nrow();
    int n1i = src_img.nrow();
    int n2 = dest_comp.nrow();
    if((n1r != n1i) || (n1i != n2))
        return;
    
    const double *src_real_pntr = &src_real(0, i);
    const double *src_img_pntr = &src_img(0, j);
    Rcomplex *dest_comp_pntr = &dest_comp(0, k);
    for(int i = 0; i < n1r; i++)
    {
        dest_comp_pntr->r = *src_real_pntr;
        dest_comp_pntr->i = *src_img_pntr;
        src_real_pntr++;
        src_img_pntr++;
        dest_comp_pntr++;
    }
}
