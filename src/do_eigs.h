#ifndef DO_EIGS_H
#define DO_EIGS_H

#include <Rcpp.h>
// R_ext/BLAS.h also includes R_ext/RS.h,
// which defines F77_CALL
#include <R_ext/BLAS.h>

extern "C" {

// ARPACK Fortran functions
void F77_CALL(dsaupd)(int *ido, char *bmat, int *n, char *which,
                      int *nev, double *tol, double *resid,
                      int *ncv, double *v, int *ldv,
                      int *iparam, int *ipntr, double *workd,
                      double *workl, int *lworkl, int *info);

void F77_CALL(dseupd)(int *rvec, char *howmny, int *select, double *d,
                      double *z, int *ldz, double *sigma, char *bmat,
                      int *n, char *which, int *nev, double *tol,
                      double *resid, int *ncv, double *v, int *ldv,
                      int *iparam, int *ipntr, double *workd, double *workl,
                      int *lworkl, int *info);

void F77_CALL(dnaupd)(int *ido, char *bmat, int *n, char *which,
                      int *nev, double *tol, double *resid,
                      int *ncv, double *v, int *ldv,
                      int *iparam, int *ipntr, double *workd,
                      double *workl, int *lworkl, int *info);

void F77_CALL(dneupd)(int *rvec, char *howmny, int *select, double *dr, double *di,
                      double *z, int *ldz, double *sigmar, double *sigmai, double *workev,
                      char *bmat, int *n, char *which, int *nev, double *tol,
                      double *resid, int *ncv, double *v, int *ldv, int *iparam,
                      int *ipntr, double *workd, double *workl, int *lworkl, int *info);

} // extern "C"

// C++ Wrapper of the functions above
inline void saupd(int& ido, char bmat, int n, char* which,
                  int nev, double& tol, double resid[],
                  int ncv, double v[], int ldv,
                  int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int& info)
{
    F77_CALL(dsaupd)(&ido, &bmat, &n, which,
                     &nev, &tol, resid,
                     &ncv, v, &ldv,
                     iparam, ipntr, workd,
                     workl, &lworkl, &info);
}

inline void seupd(bool rvec, char howmny, double d[],
                  double z[], int ldz, double sigma, char bmat,
                  int n, char *which, int nev, double tol,
                  double resid[], int ncv, double v[], int ldv,
                  int iparam[], int ipntr[], double workd[], double workl[],
                  int lworkl, int& info)
{

    int rvec_pass = (int) rvec;
    int *select_pass = new int[ncv];
    double *z_pass = (z == NULL) ? v : z;

    F77_CALL(dseupd)(&rvec_pass, &howmny, select_pass, d,
                     z_pass, &ldz, &sigma, &bmat,
                     &n, which, &nev, &tol,
                     resid, &ncv, v, &ldv,
                     iparam, ipntr, workd, workl,
                     &lworkl, &info);

    delete [] select_pass;
}

inline void naupd(int& ido, char bmat, int n, char* which,
                  int nev, double& tol, double resid[],
                  int ncv, double v[], int ldv,
                  int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int& info)
{
    F77_CALL(dnaupd)(&ido, &bmat, &n, which,
                     &nev, &tol, resid,
                     &ncv, v, &ldv,
                     iparam, ipntr, workd,
                     workl, &lworkl, &info);
}

inline void neupd(bool rvec, char howmny, double dr[], double di[],
                  double z[], int ldz, double sigmar, double sigmai, double workev[],
                  char bmat, int n, char *which, int nev, double tol,
                  double resid[], int ncv, double v[], int ldv, int iparam[],
                  int ipntr[], double workd[], double workl[], int lworkl, int& info)
{

    int rvec_pass = (int) rvec;
    int *select_pass = new int[ncv];
    double *z_pass = (z == NULL) ? v : z;

    F77_CALL(dneupd)(&rvec_pass, &howmny, select_pass, dr, di,
                     z_pass, &ldz, &sigmar, &sigmai, workev,
                     &bmat, &n, which, &nev, &tol,
                     resid, &ncv, v, &ldv, iparam,
                     ipntr, workd, workl, &lworkl, &info);

    delete [] select_pass;
}


// Prototype of function to calculate
// matrix(n, n) * vector(n) product
// A: an n*n matrix, may be sparse
// x_in: an array of length n
// y_out: result of A * x_in
// data: additional data passed to function
typedef void (*Mvfun)(SEXP mat, double *x_in, double *y_out,
                      int n, void *data);

// Common function to calculate eigen values/vectors
// mat_v_prod and data should be implemented according to
// the problem (whether matrix is dense or sparse, for example)
SEXP do_eigs_nonsym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                    SEXP params_list_r,
                    Mvfun mat_v_prod, void *data);

SEXP do_eigs_sym(SEXP A_mat_r, SEXP n_scalar_r, SEXP k_scalar_r,
                 SEXP params_list_r,
                 Mvfun mat_v_prod, void *data);


#endif // DO_EIGS_H
