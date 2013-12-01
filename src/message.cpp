#include "do_eigs.h"

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
    switch(info)
    {
        case -1:
            ::Rf_error("ARPACK/dsaupd: n must be positive");
            break;
        case -2:
            ::Rf_error("ARPACK/dsaupd: k must be positive");
            break;
        case -3:
            ::Rf_error("ARPACK/dsaupd: k < ncv <= n");
            break;
        case -4:
            ::Rf_error("ARPACK/dsaupd: maxitr must be positive");
            break;
        case -5:
            ::Rf_error("ARPACK/dsaupd: which must be one of 'LM', 'SM', 'LA', 'SA', 'BE'");
            break;
        case -6:
            ::Rf_error("ARPACK/dsaupd: error code %d", info);
            break;
        case -7:
            ::Rf_error("ARPACK/dsaupd: length of private work array WORKL is not sufficient");
            break;
        case -8:
            ::Rf_error("ARPACK/dsaupd: error return from trid. eigenvalue calculation\n"
                       "informational error from LAPACK routine dsteqr");
            break;
        case -9:
            ::Rf_error("ARPACK/dsaupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
        case -13:
            ::Rf_error("ARPACK/dsaupd: error code %d", info);
            break;
        case -9999:
            ::Rf_error("ARPACK/dsaupd: couldn't build an Arnoldi factorization");
            break;
        default:
            ::Rf_error("ARPACK/dsaupd: error code %d", info);
            break;
    }
}

// Error information
// See ARPACK/dseupd.f for details
void dseupd_error(int info)
{
    switch(info)
    {
        case -1:
            ::Rf_error("ARPACK/dseupd: n must be positive");
            break;
        case -2:
            ::Rf_error("ARPACK/dseupd: k must be positive");
            break;
        case -3:
            ::Rf_error("ARPACK/dseupd: k < ncv <= n");
            break;
        case -5:
            ::Rf_error("ARPACK/dseupd: which must be one of 'LM', 'SM', 'LA', 'SA', 'BE'");
            break;
        case -6:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
            break;
        case -7:
            ::Rf_error("ARPACK/dseupd: length of private work WORKL array is not sufficient");
            break;
        case -8:
            ::Rf_error("ARPACK/dseupd: error return from trid. eigenvalue calculation\n"
                       "informational error from LAPACK routine dsteqr");
            break;
        case -9:
            ::Rf_error("ARPACK/dseupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
            break;
        case -14:
            ::Rf_error("ARPACK/dseupd: DSAUPD did not find any eigenvalues to sufficient accuracy");
            break;
        case -15:
        case -16:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
        case -17:
            ::Rf_error("ARPACK/dseupd: DSEUPD got a different count of the number of converged Ritz values than DSAUPD got");
            break;
        default:
            ::Rf_error("ARPACK/dseupd: error code %d", info);
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
    switch(info)
    {
        case -1:
            ::Rf_error("ARPACK/dnaupd: n must be positive");
            break;
        case -2:
            ::Rf_error("ARPACK/dnaupd: k must be positive");
            break;
        case -3:
            ::Rf_error("ARPACK/dnaupd: 2 <= ncv - k <= n");
            break;
        case -4:
            ::Rf_error("ARPACK/dnaupd: maxitr must be positive");
            break;
        case -5:
            ::Rf_error("ARPACK/dnaupd: which must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'");
            break;
        case -6:
            ::Rf_error("ARPACK/dnaupd: error code %d", info);
            break;
        case -7:
            ::Rf_error("ARPACK/dnaupd: length of private work array is not sufficient");
            break;
        case -8:
            ::Rf_error("ARPACK/dnaupd: error return from LAPACK eigenvalue calculation");
            break;
        case -9:
            ::Rf_error("ARPACK/dnaupd: starting vector is zero");
            break;
        case -10:
        case -11:
        case -12:
            ::Rf_error("ARPACK/dnaupd: error code %d", info);
            break;
        case -9999:
            ::Rf_error("ARPACK/dnaupd: couldn't build an Arnoldi factorization");
            break;
        default:
            ::Rf_error("ARPACK/dnaupd: error code %d", info);
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
    switch(info)
    {
        case -1:
            ::Rf_error("ARPACK/dneupd: n must be positive");
            break;
        case -2:
            ::Rf_error("ARPACK/dneupd: k must be positive");
            break;
        case -3:
            ::Rf_error("ARPACK/dneupd: 2 <= ncv - k <= n");
            break;
        case -5:
            ::Rf_error("ARPACK/dneupd: which must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'");
            break;
        case -6:
            ::Rf_error("ARPACK/dneupd: error code %d", info);
            break;
        case -7:
            ::Rf_error("ARPACK/dneupd: length of private work WORKL array is not sufficient");
            break;
        case -8:
            ::Rf_error("ARPACK/dneupd: error return from calculation of a real Schur form");
            break;
        case -9:
            ::Rf_error("ARPACK/dneupd: error return from calculation of eigenvectors");
            break;
        case -10:
        case -11:
        case -12:
        case -13:
            ::Rf_error("ARPACK/dneupd: error code %d", info);
            break;
        case -14:
            ::Rf_error("ARPACK/dneupd: DNAUPD did not find any eigenvalues to sufficient accuracy");
        case -15:
            ::Rf_error("ARPACK/dneupd: DNEUPD got a different count of the number of converged Ritz values than DNAUPD got");
            break;
        default:
            ::Rf_error("ARPACK/dneupd: error code %d", info);
            break;
    }
}

