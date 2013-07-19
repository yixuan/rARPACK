/*
  C++ wrappers for ARPACK FORTRAN routines 
*/

#ifndef ARPACK_CPP_H
#define ARPACK_CPP_H

#include "arpack_f.h"

/* Nonsymmetric real matrix */
inline void naupp(ARint& ido, char bmat, ARint n, char* which, ARint nev,
                  double& tol, double resid[], ARint ncv, double V[],
                  ARint ldv, ARint iparam[], ARint ipntr[], double workd[],
                  double workl[], ARint lworkl, ARint& info)
{
    F77NAME(dnaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  V, &ldv, iparam, ipntr, workd, workl,
                  &lworkl, &info);
}
inline void naupp(ARint& ido, char bmat, ARint n, char* which, ARint nev,
                  float& tol, float resid[], ARint ncv, float V[],
                  ARint ldv, ARint iparam[], ARint ipntr[], float workd[],
                  float workl[], ARint lworkl, ARint& info)
{
    F77NAME(snaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  V, &ldv, iparam, ipntr, workd, workl,
                  &lworkl, &info);
}


inline void neupp(bool rvec, char HowMny, double dr[],
                  double di[], double Z[], ARint ldz, double sigmar,
                  double sigmai, double workv[], char bmat, ARint n,
                  char* which, ARint nev, double tol, double resid[],
                  ARint ncv, double V[], ARint ldv, ARint iparam[],
                  ARint ipntr[], double workd[], double workl[],
                  ARint lworkl, ARint& info)
{

  ARint      irvec;
  ARlogical* iselect;
  double*    iZ;

  irvec   = (ARint) rvec;
  iselect = new ARlogical[ncv];
  iZ = (Z == NULL) ? V : Z;

  F77NAME(dneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, workv, &bmat, &n, which, &nev, &tol,
                  resid, &ncv, V, &ldv, iparam, ipntr,
                  workd, workl, &lworkl, &info);

  delete[] iselect;
}
inline void neupp(bool rvec, char HowMny, float dr[],
                  float di[], float Z[], ARint ldz, float sigmar,
                  float sigmai, float workv[], char bmat, ARint n,
                  char* which, ARint nev, float tol, float resid[],
                  ARint ncv, float V[], ARint ldv, ARint iparam[],
                  ARint ipntr[], float workd[], float workl[],
                  ARint lworkl, ARint& info)
{

  ARint      irvec;
  ARlogical* iselect;
  float*     iZ;

  irvec   = (ARint) rvec;
  iselect = new ARlogical[ncv];
  iZ = (Z == NULL) ? V : Z;

  F77NAME(sneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, workv, &bmat, &n, which, &nev, &tol,
                  resid, &ncv, V, &ldv, iparam, ipntr,
                  workd, workl, &lworkl, &info);

  delete[] iselect;
}


/* Symmetric real matrix */
inline void saupp(ARint& ido, char bmat, ARint n, char* which, ARint nev,
                  double& tol, double resid[], ARint ncv, double V[],
                  ARint ldv, ARint iparam[], ARint ipntr[], double workd[],
                  double workl[], ARint lworkl, ARint& info)
{
    F77NAME(dsaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  V, &ldv, iparam, ipntr, workd, workl,
                  &lworkl, &info);

}
inline void saupp(ARint& ido, char bmat, ARint n, char* which, ARint nev,
                  float& tol, float resid[], ARint ncv, float V[],
                  ARint ldv, ARint iparam[], ARint ipntr[], float workd[],
                  float workl[], ARint lworkl, ARint& info)
{
    F77NAME(ssaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  V, &ldv, iparam, ipntr, workd, workl,
                  &lworkl, &info);

}


inline void seupp(bool rvec, char HowMny, double d[], double Z[],
                  ARint ldz, double sigma, char bmat, ARint n,
                  char* which, ARint nev, double tol, double resid[],
                  ARint ncv, double V[], ARint ldv, ARint iparam[],
                  ARint ipntr[], double workd[], double workl[],
                  ARint lworkl, ARint& info)
{

  ARint      irvec;
  ARlogical* iselect;
  double*    iZ;

  irvec   = (ARint) rvec;
  iselect = new ARlogical[ncv];
  iZ = (Z == NULL) ? V : Z;

  F77NAME(dseupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma, &bmat,
                  &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
                  ipntr, workd, workl, &lworkl, &info );

  delete[] iselect;

}
inline void seupp(bool rvec, char HowMny, float d[], float Z[],
                  ARint ldz, float sigma, char bmat, ARint n,
                  char* which, ARint nev, float tol, float resid[],
                  ARint ncv, float V[], ARint ldv, ARint iparam[],
                  ARint ipntr[], float workd[], float workl[],
                  ARint lworkl, ARint& info)
{

  ARint      irvec;
  ARlogical* iselect;
  float*     iZ;

  irvec   = (ARint) rvec;
  iselect = new ARlogical[ncv];
  iZ = (Z == NULL) ? V : Z;

  F77NAME(sseupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma, &bmat,
                  &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam,
                  ipntr, workd, workl, &lworkl, &info );

  delete[] iselect;

}


#endif /* ARPACK_CPP_H */
