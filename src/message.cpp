#include "do_eigs.h"
using std::string;

// Warning and error information
// See ARPACK/dsaupd.f for details
void dsaupd_warn(int info)
{
    switch(info)
    {
        case 1:
            ::Rf_warning("ARPACK/dsaupd: maximum number of iterations taken");
            break;
        case 2:
            break;
        case 3:
            ::Rf_warning("ARPACK/dsaupd: no shifts could be applied, try to increase ncv");
            break;
    }
}

void dsaupd_error(int info)
{
    char num[10];
    sprintf(num, "%d", info);
    string info_str(num);
    
    switch(info)
    {
        case -1:
            Rcpp::stop("ARPACK/dsaupd: n must be positive");
            break;
        case -2:
            Rcpp::stop("ARPACK/dsaupd: k must be positive");
            break;
        case -3:
            Rcpp::stop("ARPACK/dsaupd: k < ncv <= n");
            break;
        case -4:
            Rcpp::stop("ARPACK/dsaupd: maxitr must be positive");
            break;
        case -5:
            Rcpp::stop("ARPACK/dsaupd: which must be one of 'LM', 'SM', 'LA', 'SA', 'BE'");
            break;
        case -6:
            Rcpp::stop("ARPACK/dsaupd: error code " + info_str);
            break;
        case -7:
            Rcpp::stop("ARPACK/dsaupd: length of private work array WORKL is not sufficient");
            break;
        case -8:
            Rcpp::stop("ARPACK/dsaupd: error return from trid. eigenvalue calculation\n"
                       "informational error from LAPACK routine dsteqr");
            break;
        case -9:
            Rcpp::stop("ARPACK/dsaupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
        case -13:
            Rcpp::stop("ARPACK/dsaupd: error code " + info_str);
            break;
        case -9999:
            Rcpp::stop("ARPACK/dsaupd: couldn't build an Arnoldi factorization");
            break;
        default:
            Rcpp::stop("ARPACK/dsaupd: error code " + info_str);
            break;
    }
}

// Error information
// See ARPACK/dseupd.f for details
void dseupd_error(int info)
{
    char num[10];
    sprintf(num, "%d", info);
    string info_str(num);

    switch(info)
    {
        case -1:
            Rcpp::stop("ARPACK/dseupd: n must be positive");
            break;
        case -2:
            Rcpp::stop("ARPACK/dseupd: k must be positive");
            break;
        case -3:
            Rcpp::stop("ARPACK/dseupd: k < ncv <= n");
            break;
        case -5:
            Rcpp::stop("ARPACK/dseupd: which must be one of 'LM', 'SM', 'LA', 'SA', 'BE'");
            break;
        case -6:
            Rcpp::stop("ARPACK/dseupd: error code " + info_str);
            break;
        case -7:
            Rcpp::stop("ARPACK/dseupd: length of private work WORKL array is not sufficient");
            break;
        case -8:
            Rcpp::stop("ARPACK/dseupd: error return from trid. eigenvalue calculation\n"
                       "informational error from LAPACK routine dsteqr");
            break;
        case -9:
            Rcpp::stop("ARPACK/dseupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
            Rcpp::stop("ARPACK/dseupd: error code " + info_str);
            break;
        case -14:
            Rcpp::stop("ARPACK/dseupd: DSAUPD did not find any eigenvalues to sufficient accuracy");
            break;
        case -15:
        case -16:
            Rcpp::stop("ARPACK/dseupd: error code " + info_str);
        case -17:
            Rcpp::stop("ARPACK/dseupd: DSEUPD got a different count of the number of converged Ritz values than DSAUPD got");
            break;
        default:
            Rcpp::stop("ARPACK/dseupd: error code " + info_str);
            break;
    }
}



// Warning and error information
// See ARPACK/dnaupd.f for details
void dnaupd_warn(int info)
{
    switch(info)
    {
        case 1:
            ::Rf_warning("ARPACK/dnaupd: maximum number of iterations taken");
            break;
        case 2:
            break;
        case 3:
            ::Rf_warning("ARPACK/dnaupd: no shifts could be applied, try to increase ncv");
            break;
    }
}

void dnaupd_error(int info)
{
    char num[10];
    sprintf(num, "%d", info);
    string info_str(num);

    switch(info)
    {
        case -1:
            Rcpp::stop("ARPACK/dnaupd: n must be positive");
            break;
        case -2:
            Rcpp::stop("ARPACK/dnaupd: k must be positive");
            break;
        case -3:
            Rcpp::stop("ARPACK/dnaupd: 2 <= ncv - k <= n");
            break;
        case -4:
            Rcpp::stop("ARPACK/dnaupd: maxitr must be positive");
            break;
        case -5:
            Rcpp::stop("ARPACK/dnaupd: which must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'");
            break;
        case -6:
            Rcpp::stop("ARPACK/dnaupd: error code " + info_str);
            break;
        case -7:
            Rcpp::stop("ARPACK/dnaupd: length of private work array is not sufficient");
            break;
        case -8:
            Rcpp::stop("ARPACK/dnaupd: error return from LAPACK eigenvalue calculation");
            break;
        case -9:
            Rcpp::stop("ARPACK/dnaupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
            Rcpp::stop("ARPACK/dnaupd: error code " + info_str);
            break;
        case -9999:
            Rcpp::stop("ARPACK/dnaupd: couldn't build an Arnoldi factorization");
            break;
        default:
            Rcpp::stop("ARPACK/dnaupd: error code " + info_str);
            break;
    }
}

// Warning and error information
// See ARPACK/dneupd.f for details
void dneupd_warn(int info)
{
    switch(info)
    {
        case 1:
            ::Rf_warning("ARPACK/dneupd: the Schur form computed by LAPACK routine dlahqr"
                         "could not be reordered by LAPACK routine dtrsen");
            break;
    }
}

void dneupd_error(int info)
{
    char num[10];
    sprintf(num, "%d", info);
    string info_str(num);

    switch(info)
    {
        case -1:
            Rcpp::stop("ARPACK/dneupd: n must be positive");
            break;
        case -2:
            Rcpp::stop("ARPACK/dneupd: k must be positive");
            break;
        case -3:
            Rcpp::stop("ARPACK/dneupd: 2 <= ncv - k <= n");
            break;
        case -5:
            Rcpp::stop("ARPACK/dneupd: which must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'");
            break;
        case -6:
            Rcpp::stop("ARPACK/dneupd: error code " + info_str);
            break;
        case -7:
            Rcpp::stop("ARPACK/dneupd: length of private work WORKL array is not sufficient");
            break;
        case -8:
            Rcpp::stop("ARPACK/dneupd: error return from calculation of a real Schur form");
            break;
        case -9:
            Rcpp::stop("ARPACK/dneupd: error return from calculation of eigenvectors");
            break;
        case -10:
        case -11:
        case -12:
        case -13:
            Rcpp::stop("ARPACK/dneupd: error code " + info_str);
            break;
        case -14:
            Rcpp::stop("ARPACK/dneupd: DNAUPD did not find any eigenvalues to sufficient accuracy");
        case -15:
            Rcpp::stop("ARPACK/dneupd: DNEUPD got a different count of the number of converged Ritz values than DNAUPD got");
            break;
        default:
            Rcpp::stop("ARPACK/dneupd: error code " + info_str);
            break;
    }
}

