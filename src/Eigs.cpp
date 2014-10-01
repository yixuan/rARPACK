#include "Eigs.h"

Eigs::Eigs(int n_, int nev_, int ncv_, MatOp *op_,
           const string & which_, int workmode_,
           char bmat_, double tol_, int maxitr_)
{
    n = n_;
    nev = nev_;
    ncv = ncv_;
    op = op_;
    which = which_;
    workmode = workmode_;
    bmat = bmat_;
    tol = tol_;
    maxitr = maxitr_;

    ido = 0;
    info = 0;
    ierr = 0;
    updatecount = 0;

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

void Eigs::update()
{
    initResid();

    while (ido != 99)
    {
        aupd();
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

void Eigs::checkUpdateError()
{
    // Ensure that update() is called at least once
    if (updatecount < 1)
        Rcpp::stop("need to call update() first");

    // info > 0 means warning, < 0 means error
    if (info > 0)  warning(1, info);
    if (info < 0)  error(1, info);
}

Eigs::~Eigs()
{
    delete [] workd;
    delete [] resid;
}

