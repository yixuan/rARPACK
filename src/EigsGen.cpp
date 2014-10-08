#include "EigsGen.h"
#include <RcppEigen.h>

using std::string;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::RealSchur;
using Eigen::EigenSolver;

typedef Eigen::Map<VectorXd> MapVec;
typedef Eigen::Map<MatrixXd> MapMat;

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
    wl = new double[lworkl];
    vm = new double[n * ncv];
}


EigsGen::~EigsGen()
{
    delete [] workv;
    delete [] workl;
    delete [] wl;
    delete [] vm;
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

void EigsGen::aupd()
{
    naupd(ido, bmat, n, which.c_str(),
          nev, tol, resid,
          ncv, eigV.begin(), n,
          iparam, ipntr, workd,
          workl, lworkl, info);
}

void EigsGen::eupd()
{
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Used to store results, will use V instead.
    double *Z = eigV.begin();
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    
    std::copy(workl, workl + lworkl, wl);
    std::copy(eigV.begin(), eigV.end(), vm);
    
    // Use neupd() to retrieve results
    neupd(retvec, howmny, eigdr.begin(), eigdi.begin(),
          Z, ldz, op->getsigmar(), op->getsigmai(), workv,
          bmat, n, which.c_str(), nev, tol,
          resid, ncv, eigV.begin(), n,
          iparam, ipntr, workd, workl,
          lworkl, ierr);
}

std::complex<double> EigsGen::eigenvalue2x2(const double &a,
    const double &b, const double &c, const double &d)
{
    double real = (a + d) * 0.5;
    double imag = 0.5 * sqrt(4 * (a * d - b * c) - (a + d) * (a + d));
    return std::complex<double>(real, imag);
}

void EigsGen::eigenvalueSchur(const MatrixXd &Rm, VectorXcd &result)
{
    int m = Rm.cols();
    if(result.size() != m)
        result.resize(m);
    for(int i = 0; i < m; i++)
    {
        if(i == m - 1 || fabs(Rm(i + 1, i)) < 1e-16)
        {
            result[i] = std::complex<double>(Rm(i, i), 0);
        } else {
            result[i] = eigenvalue2x2(Rm(i, i), Rm(i, i + 1),
                                      Rm(i + 1, i), Rm(i + 1, i + 1));
            i++;
            result[i] = conj(result[i - 1]);
        }
    }
}

void EigsGen::findMatchedIndex(const Eigen::VectorXcd &target,
                               const Eigen::VectorXcd &collection,
                               Eigen::VectorXi &result)
{
    int nfound = 0;
    int maxn = target.size();
    if(result.size() < maxn)
        result.resize(maxn);
    for(int i = 0; i < collection.size(); i++)
    {
        int j;
        for(j = 0; j < maxn; j++)
        {
            if(abs(collection[i] - target[j]) < 1e-10)
                break;
        }
        if(j < maxn)
        {
            result[nfound] = i;
            nfound++;
            if(collection[i].imag() != 0)
            {
                i++;
                result[nfound] = i;
                nfound++;
            }
        }
        if(nfound >= maxn)  break;
    }
    if(result.size() > nfound)
        result.conservativeResize(nfound);
}

Rcpp::List EigsGen::extract2()
{
    Rcpp::NumericMatrix H(ncv, ncv);
    Rcpp::NumericMatrix V(n, ncv);
    Rcpp::NumericVector real(ncv);
    Rcpp::NumericVector imag(ncv);
    int nconv = iparam[5 - 1];
    int niter = iparam[9 - 1];
    
    std::copy(wl, wl + ncv * ncv, H.begin());
    std::copy(vm, vm + n * ncv, V.begin());
    std::copy(wl + ncv * ncv, wl + ncv * ncv + ncv, real.begin());
    std::copy(wl + ncv * ncv + ncv, wl + ncv * ncv + 2 * ncv, imag.begin());

    MapMat Hm(wl, ncv, ncv);
    MapMat Vm(vm, n, ncv);
    RealSchur<MatrixXd> schur(Hm);
    MatrixXd Qm = schur.matrixU();
    MatrixXd Rm = schur.matrixT();
    VectorXcd evalsRm(ncv);
    VectorXi selectInd(nconv);
    
    VectorXcd evalsConverged(nconv);
    evalsConverged.real() = MapVec(wl + ncv * ncv, nconv);
    evalsConverged.imag() = MapVec(wl + ncv * ncv + ncv, nconv);
    
    eigenvalueSchur(Rm, evalsRm);
    findMatchedIndex(evalsConverged.head(nconv), evalsRm, selectInd);
    
    int nfound = selectInd.size();
    if(nfound < 1)  return R_NilValue;
    
    // Shrink Qm and Rm to the dimension given by the largest value
    // in selectInd. Since selectInd is strictly increasing,
    // we can just use its last value.
    int lastInd = selectInd[selectInd.size() - 1];
    Qm.conservativeResize(Eigen::NoChange, lastInd + 1);
    Rm.conservativeResize(lastInd + 1, lastInd + 1);
    
    // Eigen decomposition of Rm
    EigenSolver<MatrixXd> es(Rm);
    evalsRm = es.eigenvalues();
    MatrixXcd eigenvectors = Vm * (Qm * es.eigenvectors());
    for(int i = 0; i < nfound; i++)
    {
        // Since selectInd[i] >= i for all i, it is safe to
        // overwrite the elements and columns.
        evalsRm[i] = evalsRm[selectInd[i]];
        eigenvectors.col(i) = eigenvectors.col(selectInd[i]);
    }
    evalsRm.conservativeResize(nfound);
    eigenvectors.conservativeResize(Eigen::NoChange, nfound);
    
    return Rcpp::List::create(Rcpp::Named("values") = evalsRm,
                              Rcpp::Named("vectors") = eigenvectors,
                              Rcpp::Named("V") = V,
                              Rcpp::Named("H") = H,
                              Rcpp::Named("Q") = Qm,
                              Rcpp::Named("R") = Rm,
                              Rcpp::Named("ind") = selectInd,
                              Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                              Rcpp::Named("niter") = Rcpp::wrap(niter));
}

Rcpp::List EigsGen::extract()
{
    // Result list
    Rcpp::List ret;
    // Obtain 'nconv' converged eigenvalues
    int nconv = iparam[5 - 1];
    // 'niter' number of iterations
    int niter = iparam[9 - 1];

    /*********************************************/
    //
    // Case 1: If there is no converged eigenvalue 
    //
    /*********************************************/
    if(nconv <= 0)
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
    if(nconv < nev)
        ::Rf_warning("only %d eigenvalues converged, less than k", nconv);
    // Sometimes there are nconv = nev + 1 converged eigenvalues,
    // mainly due to pairs of complex eigenvalues.
    // We will truncate at nev
    int truenconv = nconv > nev ? nev : nconv;

    // equivalent R code: if(all(abs(dimag[1:nconv] < 1e-17)))
    if(Rcpp::is_true(Rcpp::all(Rcpp::abs(eigdi) < 1e-17)))
    {
        // v.erase(start, end) removes v[start <= i < end]
        eigdr.erase(truenconv, eigdr.length());
        // Make sure eigenvalues are in decreasing order.
        Rcpp::IntegerVector order = sort_with_order(eigdr);
        
        if(!retvec)
        {
            ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                     Rcpp::Named("vectors") = R_NilValue,
                                     Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                                     Rcpp::Named("niter") = Rcpp::wrap(niter));
            return ret;
        }
        
        // The matrix to be returned
        Rcpp::NumericMatrix retV(n, truenconv);
        for(int i = 0; i < truenconv; i++)
        {
            copy_column(eigV, order[i], retV, i);
        }

        ret = Rcpp::List::create(Rcpp::Named("values") = eigdr,
                                 Rcpp::Named("vectors") = retV,
                                 Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));

        return ret;
    }

    /*************************************/
    //
    // Case 3: Complex eigenvalues
    //
    /*************************************/
    Rcpp::ComplexVector cmpeigd(nconv);
    for(int i = 0; i < nconv; i++)
    {
        cmpeigd[i].r = eigdr[i];
        cmpeigd[i].i = eigdi[i];
    }
    
    if(!retvec)
    {
        // To make sure the eigenvalues are in proper order
        sort_with_order(cmpeigd);
        if(nconv > truenconv)  cmpeigd.erase(truenconv, nconv);
        ret = Rcpp::List::create(Rcpp::Named("values") = cmpeigd,
                                 Rcpp::Named("vectors") = R_NilValue,
                                 Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                                 Rcpp::Named("niter") = Rcpp::wrap(niter));
        return ret;
    }
    
    Rcpp::ComplexMatrix retV(n, truenconv);
    Rcpp::ComplexVector retd = Rcpp::clone(cmpeigd);
    bool *img_zero = new bool[nconv + 1];  // add 1 to be safe
    for(int i = 0; i < nconv; i++)
    {
        if(fabs(eigdi[i]) < 1e-16)
        {
            img_zero[i] = true;
        } else {
            img_zero[i] = false;
            i++;
            img_zero[i] = false;
            retd[i].r = 0;
            retd[i].i = 0;
        }
    }
    Rcpp::IntegerVector order = sort_with_order(retd);
    int *order_pntr = order.begin();
    for(int i = 0; i < truenconv; i++, order_pntr++)
    {
        int ind = *order_pntr;
        if(img_zero[ind])
        {
            retd[i] = cmpeigd[ind];
            copy_column(eigV, ind, eigV, ind, 0, retV, i);
        } else {
            retd[i] = cmpeigd[ind];
            copy_column(eigV, ind, eigV, ind + 1, 1, retV, i);
            i++;
            if(i < truenconv)
            {
                retd[i].r = cmpeigd[ind].r;
                retd[i].i = -cmpeigd[ind].i;
                copy_column(eigV, ind, eigV, ind + 1, -1, retV, i);
            }
        }
    }
    delete [] img_zero;
    if(nconv > truenconv)  retd.erase(truenconv, nconv);
    
    ret = Rcpp::List::create(Rcpp::Named("values") = retd,
                             Rcpp::Named("vectors") = retV,
                             Rcpp::Named("nconv") = Rcpp::wrap(truenconv),
                             Rcpp::Named("niter") = Rcpp::wrap(niter),
                             extract2());

    return ret;
}
