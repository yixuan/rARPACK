#include "EigsSym.h"

EigsSym::EigsSym(int n_, int nev_, int ncv_, MatOp *op_,
                 const string & which_, int workmode_,
                 char bmat_, double tol_, int maxitr_) :
    Eigs(n_, nev_, ncv_, op_, which_, workmode_,
        bmat_, tol_, maxitr_),
    eigV(n, ncv), eigd(nev)
{
    lworkl = ncv * (ncv + 8);
    workl = new double[lworkl]();
    updatecount = 0;
}


EigsSym::~EigsSym()
{
    delete [] workl;
    delete [] workd;
    delete [] resid;
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

void EigsSym::update()
{
    initResid();

    while (ido != 99)
    {
        saupd(ido, bmat, n, which.c_str(),
              nev, tol, resid,
              ncv, eigV.begin(), n,
              iparam, ipntr, workd,
              workl, lworkl, info);
        switch(ido)
        {
            case -1:
            case 1:
                // Shift-and-invert
                if (workmode == 3)
                    op->shiftSolve(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
                else
                    op->prod(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
                break;
            default:
                break;
        }
        updatecount++;
    }
}

void EigsSym::checkUpdateError()
{
    // Ensure that update() is called at least once
    if (updatecount < 1)
        Rcpp::stop("need to call Update() first");

    // info > 0 means warning, < 0 means error
    if (info > 0)  warning(1, info);
    if (info < 0)  error(1, info);
}

Rcpp::List EigsSym::extract(bool rvec)
{
    checkUpdateError();

    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Used to store results, will use V instead.
    double *Z = eigV.begin();
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;

    // Number of converged eigenvalues
    int nconv = 0;
    // Number of iterations
    int niter = 0;
    // Result list
    Rcpp::List ret;

    // Use seupd() to retrieve results
    seupd(rvec, howmny, eigd.begin(),
          Z, ldz, op->getsigmar(), bmat,
          n, which.c_str(), nev, tol,
          resid, ncv, eigV.begin(), n,
          iparam, ipntr, workd, workl,
          lworkl, ierr);

    // Obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    // 'niter' number of iterations
    niter = iparam[9 - 1];

    // ierr < 0 means error
    if (ierr < 0)  error(2, ierr);

    if (nconv <= 0)
    {
        ::Rf_warning("no converged eigenvalues found");
        ret = Rcpp::List::create(Rcpp::Named("values") = R_NilValue,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
    } else {
        if (nconv < nev)
            ::Rf_warning("only %d eigenvalues converged, less than k", nconv);
        // v.erase(start, end) removes v[start <= i < end]
        eigd.erase(nconv, eigd.length());
        // ARPACK gives eigenvalues in increasing order.
        // We need decreasing one.
        std::reverse(eigd.begin(), eigd.end());
        if (rvec)
        {
            // Also change the order of columns of v_ret
            for(int i = 0; i < nconv / 2; i++)
            {
                std::swap_ranges(&eigV(0, i), &eigV(0, i + 1),
                                 &eigV(0, nconv - i - 1));
            }
            Rcpp::Range range = Rcpp::Range(0, nconv - 1);
            ret = Rcpp::List::create(Rcpp::Named("values") = eigd,
                                     Rcpp::Named("vectors") = eigV(Rcpp::_, range),
                                     Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
        } else {
            ret = Rcpp::List::create(Rcpp::Named("values") = eigd,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
        }
    }

    return ret;
}
