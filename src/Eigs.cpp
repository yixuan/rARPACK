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



// Sort the array and return the order
typedef std::pair<double, int> ValInd;
inline bool compare_ascend(const ValInd &l, const ValInd &r)
{
    return l.first < r.first;
}
inline bool compare_descend(const ValInd &l, const ValInd &r)
{
    return l.first > r.first;
}
Rcpp::IntegerVector sort_with_order(Rcpp::NumericVector &array,
                                    bool descend)
{
    int len = array.length();
    Rcpp::IntegerVector order(len);
    double *valptr = array.begin();
    int *indptr = order.begin();
    
    std::vector<ValInd> v(len);
    for(int i = 0; i < len; i++)
    {
        v[i].first = valptr[i];
        v[i].second = i;
    }
    if(descend)
        std::sort(v.begin(), v.end(), compare_descend);
    else
        std::sort(v.begin(), v.end(), compare_ascend);
    
    for(int i = 0; i < len; i++)
    {
        valptr[i] = v[i].first;
        indptr[i] = v[i].second;
    }
    
    return order;
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
