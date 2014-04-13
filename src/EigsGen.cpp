#include "EigsGen.h"

EigsGen::EigsGen(int n_, int nev_, int ncv_, MatOp *op_,
                 const string & which_, int workmode_,
                 char bmat_, double tol_, int maxitr_) :
    Eigs(n_, nev_, ncv_, op_, which_, workmode_,
         bmat_, tol_, maxitr_),
    eigV(n, ncv), eigdr(nev + 1), eigdi(nev + 1)
{
    lworkl = 3 * ncv * ncv + 6 * ncv;
    workl = new double[lworkl]();
    workv = new double[3 * ncv]();
    updatecount = 0;
}


EigsGen::~EigsGen()
{
    delete [] workv;
    delete [] workl;
    delete [] workd;
    delete [] resid;
}

void EigsGen::error(int stage, int errorcode)
{
    if (stage == 1) // dnaupd
    {
        dnaupd_error(errorcode);
    } else { // dneupd
        dneupd_error(errorcode);
    }
}

void EigsGen::warning(int stage, int errorcode)
{
    if (stage == 1) // dnaupd
    {
        dnaupd_warn(errorcode);
    } else { // dneupd
        dneupd_warn(errorcode);
    }
}

void EigsGen::update()
{
    initResid();

    while (ido != 99)
    {
        naupd(ido, bmat, n, which.c_str(),
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

void EigsGen::checkUpdateError()
{
    // Ensure that update() is called at least once
    if (updatecount < 1)
        Rcpp::stop("need to call Update() first");

    // info > 0 means warning, < 0 means error
    if (info > 0)  warning(1, info);
    if (info < 0)  error(1, info);
}

Rcpp::List EigsGen::extract(bool rvec)
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
    neupd(rvec, howmny, eigdr.begin(), eigdi.begin(),
          Z, ldz, op->getsigmar(), op->getsigmai(), workv,
          bmat, n, which.c_str(), nev, tol,
          resid, ncv, eigV.begin(), n,
          iparam, ipntr, workd, workl,
          lworkl, ierr);

    // Obtain 'nconv' converged eigenvalues
    nconv = iparam[5 - 1];
    // 'niter' number of iterations
    niter = iparam[9 - 1];

    // ierr > 0 means warning, < 0 means error
    if (ierr > 0) warning(2, ierr);
    if (ierr < 0) error(2, ierr);

    /*********************************************/
    //
    // Case 1: If there is no converged eigenvalue 
    //
    /*********************************************/
    if (nconv <= 0)
    {
        ::Rf_warning("no converged eigenvalues found");
        ret = Rcpp::List::create(Rcpp::Named("values") = R_NilValue,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
        return ret;
    }

    /*************************************/
    //
    // Case 2: If all eigenvalues are real
    //
    /*************************************/
    if (nconv < nev)
        ::Rf_warning("only %d eigenvalues converged, less than k", nconv);
    // Sometimes there are nconv = nev + 1 converged eigenvalues,
    // mainly due to pairs of complex eigenvalues.
    // We will truncate at nev
    if (nconv > nev)  nconv = nev;

    // equivalent R code: if (all(abs(dimag[1:nconv] < 1e-17)))
    if (Rcpp::is_true(Rcpp::all(Rcpp::abs(eigdi) < 1e-17)))
    {
        // v.erase(start, end) removes v[start <= i < end]
        eigdr.erase(nconv, eigdr.length());
        if (rvec)
        {
            Rcpp::Range range = Rcpp::Range(0, nconv - 1);
            ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                     Rcpp::Named("vectors") = eigV(Rcpp::_, range),
                                     Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
        } else {
            ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
        }
        return ret;
    }

    /*************************************/
    //
    // Case 3: Complex eigenvalues
    //
    /*************************************/
    Rcpp::ComplexVector cmpeigd(nconv);
    for (int i = 0; i < nconv; i++)
    {
        cmpeigd[i].r = eigdr[i];
        cmpeigd[i].i = eigdi[i];
    }
    if (!rvec)
    {
        // It seems that when rvec is false, ARPACK will return
        // the eigenvalues in ascending order. So here we reverse
        // the vector
        std::reverse(cmpeigd.begin(), cmpeigd.end());
        ret = Rcpp::List::create(Rcpp::Named("values") = cmpeigd,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
    } else {
        Rcpp::ComplexMatrix cmpeigV(n, nconv);
        // Obtain the real and imaginary part of the eigenvectors
        bool first = true;
        for (int i = 0; i < nconv; i++)
        {
            if (fabs(eigdi[i]) > 1e-30)
            {
                if (first)
                {
                    for (int j = 0; j < n; j++)
                    {
                        cmpeigV(j, i).r = eigV(j, i);
                        cmpeigV(j, i).i = eigV(j, i + 1);
                        if (i + 1 < nconv)
                            cmpeigV(j, i + 1).r = eigV(j, i);
                        if (i + 1 < nconv)
                            cmpeigV(j, i + 1).i = -eigV(j, i + 1);
                    }
                    first = false;
                } else {
                    first = true;
                }
            } else {
                for (int j = 0; j < n; j++)
                {
                    cmpeigV(j, i).r = eigV(j, i);
                    cmpeigV(j, i).i = 0;
                }
                first = true;
            }
        }
        ret = Rcpp::List::create(Rcpp::Named("values") = cmpeigd,
                                 Rcpp::Named("vectors") = cmpeigV,
                                 Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
    }

    return ret;
}
