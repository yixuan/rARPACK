#include "Eigs.h"

using std::string;

Eigs::Eigs(int n_, int nev_, int ncv_, MatOp *op_,
           const string & which_, int workmode_,
           char bmat_, double tol_, int maxitr_) :
    bmat(bmat_), n(n_), which(which_), nev(nev_), tol(tol_),
    ncv(ncv_), maxitr(maxitr_),
    workmode(workmode_), op(op_),
    ido(0), info(0), ierr(0), retvec(false)
{
    for(int i = 0; i < 11; i++)
        iparam[i] = 0;

    iparam[1 - 1] = 1;
    iparam[3 - 1] = maxitr;
    iparam[4 - 1] = 1;
    iparam[7 - 1] = workmode;

    for(int i = 0; i < 14; i++)
        ipntr[i] = 0;

    resid = new double[n]();
    workd = new double[3 * n]();
    // lworkl and workl are not initialized here since their
    // sizes are determined according to the problem.
    // Classes derived from Eigs should do this.
}

void Eigs::initResid()
{
    // Create initial residual vector
    // info = 1 means using the residual vector we provide
    info = 1;
    #include "rands.h"
    if(n <= rands_len)
    {
        op->prod(rands, resid);
    } else {
        double *initcoef = new double[n];
        double *coef_pntr = initcoef;
        for(int i = 0; i < n / rands_len; i++, coef_pntr += rands_len)
        {
            std::copy(rands, rands + rands_len, coef_pntr);
        }
        std::copy(rands, rands + n % rands_len, coef_pntr);
        // resid = A * initcoef
        op->prod(initcoef, resid);
        delete [] initcoef;
    }
}

void Eigs::matOp(double *x_in, double *y_out)
{
    if (workmode == 3)  // Shift-and-invert mode
        op->shiftSolve(x_in, y_out);
    else                // Regular mode
        op->prod(x_in, y_out);
}

void Eigs::matOp()
{
    matOp(&workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
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
