#ifndef EIGS_H
#define EIGS_H

#include <Rcpp.h>
#include "MatOp.h"


// Need to be implemented:
//     error()
//     warning()
//     aupd()
//     eupd()
//     extract()
//
class Eigs
{
protected:
    // Parameters to be set
    //
    // 'I' means standard eigenvalue problem,
    //     A * x = lambda * x
    // 'G' for generalized eigenvalue problem,
    //     A * x = lambda * B * x
    char bmat;
    // Dimension of A (n x n)
    int n;
    // Specify selection criteria
    // "LM": largest magnitude
    // "SM": smallest magnitude
    // "LR", "LI": largest real/imaginary part
    // "SR", "SI": smallest real/imaginary part
    // "LA": largest (algebraic) eigenvalues
    // "SA": smallest (algebraic) eigenvalues
    // "BE": half of each end
    std::string which;
    // Number of eigenvalues requested
    int nev;
    // Precision
    double tol;
    // Related to the algorithm, large ncv results in
    // faster convergence, but with greater memory use
    int ncv;
    // Maximum number of iterations
    int maxitr;
    // Workmode
    // workmode = 1: regular mode, we only need matrix product
    // workmode = 3: shift-and-invert mode,
    //               need to solve linear equation
    int workmode;

    // Matrix operation object
    MatOp *op;

    // Variables for computation in ARPACK
    //
    // Convergence flag
    int ido;
    // Warning and error flag in the first stage
    int info;
    // Warning and error flag in the second stage
    int ierr;
    // Control parameters
    int iparam[11];
    // Some pointers
    int ipntr[14];
    // Residual vector
    double *resid;
    // Working space and corresponding dimensions
    double *workd;
    int lworkl;
    double *workl;
    // Whether to return eigenvector
    bool retvec;

    // stage = 1: _aupd
    // stage = 2: _eupd
    virtual void error(int stage, int errorcode) = 0;
    virtual void warning(int stage, int errorcode) = 0;
    
    // Mimic C RNG
    // http://stackoverflow.com/questions/4768180/rand-implementation
    static unsigned long int seed_next;
    void srand(unsigned int seed);
    unsigned int rand(); // RAND_MAX assumed to be 32767
    
    // Generate initial residual vector
    void initResid();
    
    // Wrapper of _aupd and _eupd
    virtual void aupd() = 0;
    virtual void eupd() = 0;
    
    // matrix operation
    void matOp();
public:
    // Constructor
    Eigs(int n_, int nev_, int ncv_, MatOp *op_,
         const std::string & which_ = "LM", int workmode_ = 1,
         char bmat_ = 'I', double tol_ = 1e-10, int maxitr_ = 1000);
    // ARPACK computing
    void compute(bool rvec = true);
    // Extract results
    virtual Rcpp::List extract() = 0;
    // Destructor
    virtual ~Eigs();
};



// Helper functions to assist sorting eigenvalues
template<int RTYPE>
struct sort_order
{
    typedef typename Rcpp::traits::storage_type<RTYPE>::type type;
    typedef std::pair<type, int> ValInd;
};

// desc = 0 means ascending
// desc = 1 means descending
enum SORTORDER { ASCEND = 0, DESCEND };

template<int RTYPE, int desc>
inline bool compare_val(const typename sort_order<RTYPE>::ValInd &l,
                        const typename sort_order<RTYPE>::ValInd &r)
{
    if(desc)
        return l.first > r.first;
    else
        return l.first < r.first;
}
template<>
inline bool compare_val<CPLXSXP, ASCEND>(const sort_order<CPLXSXP>::ValInd &l,
                                         const sort_order<CPLXSXP>::ValInd &r)
{
    return l.first.r * l.first.r + l.first.i * l.first.i <
           r.first.r * r.first.r + r.first.i * r.first.i;
}
template<>
inline bool compare_val<CPLXSXP, DESCEND>(const sort_order<CPLXSXP>::ValInd &l,
                                          const sort_order<CPLXSXP>::ValInd &r)
{
    return l.first.r * l.first.r + l.first.i * l.first.i >
           r.first.r * r.first.r + r.first.i * r.first.i;
}

// Sort the array and return the order
template<int RTYPE>
inline Rcpp::IntegerVector sort_with_order(Rcpp::Vector<RTYPE> &array,
                                           bool descend = true)
{
    int len = array.length();
    Rcpp::IntegerVector order(len);
    typename sort_order<RTYPE>::type *valptr = array.begin();
    int *indptr = order.begin();
    
    std::vector<typename sort_order<RTYPE>::ValInd> v(len);
    for(int i = 0; i < len; i++)
    {
        v[i].first = valptr[i];
        v[i].second = i;
    }
    if(descend)
        std::sort(v.begin(), v.end(), compare_val<RTYPE, DESCEND>);
    else
        std::sort(v.begin(), v.end(), compare_val<RTYPE, ASCEND>);
    
    for(int i = 0; i < len; i++)
    {
        valptr[i] = v[i].first;
        indptr[i] = v[i].second;
    }
    
    return order;
}

void copy_column(const Rcpp::NumericMatrix &source, int i,
                 Rcpp::NumericMatrix &dest, int j);
void copy_column(const Rcpp::NumericMatrix &src_real, int i,
                 const Rcpp::NumericMatrix &src_img, int j, int sign,
                 Rcpp::ComplexMatrix &dest_comp, int k);

#endif // EIGS_H

