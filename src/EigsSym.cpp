#include "EigsSym.h"

using std::string;

EigsSym::EigsSym(int n_, int nev_, int ncv_, MatOp *op_,
                 const string & which_, int workmode_,
                 char bmat_, double tol_, int maxitr_) :
    Eigs(n_, nev_, ncv_, op_, which_, workmode_,
        bmat_, tol_, maxitr_),
    eigV(n, ncv), eigd(nev)
{
    lworkl = ncv * (ncv + 8);
    workl = new double[lworkl]();
}


EigsSym::~EigsSym()
{
    delete [] workl;
}

void EigsSym::error(int stage, int errorcode)
{
    if (stage == 1) // dsaupd
    {
        dsaupd_error(errorcode);
    } else { // dseupd
        dseupd_error(errorcode);
    }
}

void EigsSym::warning(int stage, int errorcode)
{
    if (stage == 1) // dsaupd
    {
        dsaupd_warn(errorcode);
    } else { // dseupd
        // no warning for dsuepd
    }
}

void EigsSym::aupd()
{
    saupd(ido, bmat, n, which.c_str(),
          nev, tol, resid,
          ncv, eigV.begin(), n,
          iparam, ipntr, workd,
          workl, lworkl, info);
}

void EigsSym::eupd()
{
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Used to store results, will use V instead.
    double *Z = eigV.begin();
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    
    // Use seupd() to retrieve results
    seupd(retvec, howmny, eigd.begin(),
          Z, ldz, op->getsigmar(), bmat,
          n, which.c_str(), nev, tol,
          resid, ncv, eigV.begin(), n,
          iparam, ipntr, workd, workl,
          lworkl, ierr);
}



typedef std::pair<double, int> ValInd;
// SORTORDER = 0 means ascending
// SORTORDER = 1 means descending
enum { ASCEND = 0, DESCEND };

template<int SORTORDER>
bool compare_val(const ValInd &l, const ValInd &r)
{
    if(SORTORDER == DESCEND)
        return l.first > r.first;
    else
        return l.first < r.first;
}

// Sort the array and return the order
template<int SORTORDER>
Rcpp::IntegerVector sort_with_order(Rcpp::NumericVector &array)
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
    std::sort(v.begin(), v.end(), compare_val<SORTORDER>);
    
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



Rcpp::List EigsSym::extract()
{
    // Result list
    Rcpp::List ret;
    // Obtain 'nconv' converged eigenvalues
    int nconv = iparam[5 - 1];
    // 'niter' number of iterations
    int niter = iparam[9 - 1];

    if(nconv <= 0)
    {
        ::Rf_warning("no converged eigenvalues found");
        return returnResult(R_NilValue, R_NilValue,
                            Rcpp::wrap(nconv), Rcpp::wrap(niter));
    }

    if(nconv < nev)
        ::Rf_warning("only %d eigenvalues converged, less than k", nconv);
    
    // v.erase(start, end) removes v[start <= i < end]
    eigd.erase(nconv, eigd.length());
    // ARPACK gives eigenvalues in increasing order.
    // We need decreasing one.
    Rcpp::IntegerVector order = sort_with_order<DESCEND>(eigd);
    
    if(!retvec)
    {
        return returnResult(eigd, R_NilValue, Rcpp::wrap(nconv),
                            Rcpp::wrap(niter));
    }
    
    // The matrix to be returned
    Rcpp::NumericMatrix retV(n, nconv);
    for(int i = 0; i < nconv; i++)
    {
        copy_column(eigV, order[i], retV, i);
    }

    return returnResult(eigd, retV, Rcpp::wrap(nconv),
                        Rcpp::wrap(niter));
}
