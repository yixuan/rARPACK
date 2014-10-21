#ifndef EIGSSYM_H
#define EIGSSYM_H

#include "Eigs.h"
#include "ARPACK.h"


// Implemented:
//     error()
//     warning()
//     aupd()
//     eupd()
//     extract()
//
class EigsSym: public Eigs
{
protected:
    virtual void error(int stage, int errorcode);
    virtual void warning(int stage, int errorcode);

    // Store final results of eigenvectors
    // Conceptually it is an n * ncv matrix
    Rcpp::NumericMatrix eigV;
    // Store final results of eigenvalues
    // Conceptually it is a vector of length nev
    Rcpp::NumericVector eigd;

    virtual void aupd();
    virtual void eupd();
public:
    EigsSym(int n_, int nev_, int ncv_, MatOp *op_,
            const std::string & which_ = "LM", int workmode_ = 1, 
            char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    virtual Rcpp::List extract();
    virtual ~EigsSym();
};



// Helper functions to sort eigenvalues

typedef std::pair<double, int> ValInd;
// SORTORDER = 0 means ascending
// SORTORDER = 1 means descending
enum { ASCEND = 0, DESCEND };

template<int SORTORDER>
inline bool compare_val(const ValInd &l, const ValInd &r)
{
    if(SORTORDER == DESCEND)
        return l.first > r.first;
    else
        return l.first < r.first;
}

// Sort the array and return the order
template<int SORTORDER>
inline Rcpp::IntegerVector sort_with_order(Rcpp::NumericVector &array)
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
inline void copy_column(const Rcpp::NumericMatrix &source, int i,
                        Rcpp::NumericMatrix &dest, int j)
{
    int n1 = source.nrow();
    int n2 = dest.nrow();
    if(n1 != n2)  return;
    
    std::copy(&source(0, i), &source(0, i) + n1, &dest(0, j));
}

#endif // EIGSSYM_H

