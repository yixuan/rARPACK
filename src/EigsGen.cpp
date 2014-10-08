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
         bmat_, tol_, maxitr_)
{
    V = new double[n * ncv];
    dr = new double[nev + 1];
    di = new double[nev + 1];
    
    lworkl = 3 * ncv * ncv + 6 * ncv;
    workl = new double[lworkl]();
    workv = new double[3 * ncv]();
}


EigsGen::~EigsGen()
{
    delete [] workv;
    delete [] workl;
    delete [] di;
    delete [] dr;
    delete [] V;
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
          ncv, V, n,
          iparam, ipntr, workd,
          workl, lworkl, info);
}

void EigsGen::eupd()
{
    // neupd() sometimes has bugs (e.g., calculating wrong
    // Hessenburg matrix)
    // We try to implement neupd() by our own. This is done in
    // extract().
    
    // Below is the original code
    
    /*
    // 'A' means to calculate Ritz vectors
    // 'P' to calculate Schur vectors
    char howmny = 'A';
    // Used to store results, will use V instead.
    double *Z = eigV.begin();
    // Leading dimension of Z, required by FORTRAN
    int ldz = n;
    
    // Use neupd() to retrieve results
    neupd(retvec, howmny, eigdr.begin(), eigdi.begin(),
          Z, ldz, op->getsigmar(), op->getsigmai(), workv,
          bmat, n, which.c_str(), nev, tol,
          resid, ncv, eigV.begin(), n,
          iparam, ipntr, workd, workl,
          lworkl, ierr);
    */
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

void EigsGen::recomputeH()
{
    MapMat Hm(workl, ncv, ncv);
    MapMat Vm(V, n, ncv);
    MatrixXd AV(n, ncv);
    for(int i = 0; i < ncv; i++)
    {
        matOp(&Vm(0, i), &AV(0, i));
    }
    Hm = Vm.topRows(ncv).householderQr().solve(AV.topRows(ncv));
}

void EigsGen::transformEigenvalues(VectorXcd &evals)
{
    // Only transform eigenvalues when using real sigma
    if(workmode == 3 && fabs(op->getsigmai()) < 1e-16)
        evals = 1.0 / evals.array() + op->getsigmar();
}

Rcpp::List EigsGen::extract()
{
    recomputeH();

    int nconv = iparam[5 - 1];
    int niter = iparam[9 - 1];

    MapMat Hm(workl, ncv, ncv);
    MapMat Vm(V, n, ncv);
    RealSchur<MatrixXd> schur(Hm);
    MatrixXd Qm = schur.matrixU();
    MatrixXd Rm = schur.matrixT();
    VectorXcd evalsRm(ncv);
    VectorXi selectInd(nconv);
    
    VectorXcd evalsConverged(nconv);
    evalsConverged.real() = MapVec(workl + ncv * ncv, nconv);
    evalsConverged.imag() = MapVec(workl + ncv * ncv + ncv, nconv);
    
    eigenvalueSchur(Rm, evalsRm);
    findMatchedIndex(evalsConverged.head(nconv), evalsRm, selectInd);
    
    //Rcpp::Rcout << evalsRm << "\n\n";
    //Rcpp::Rcout << evalsConverged << "\n\n";
    //::Rf_PrintValue(Rcpp::wrap(MapMat(wl + ipntr[4] - 1, ncv, ncv)));
    
    int nfound = selectInd.size();
    if(nfound < 1)
    {
        ::Rf_warning("no converged eigenvalues found");
        return Rcpp::List::create(
            Rcpp::Named("values") = R_NilValue,
            Rcpp::Named("vectors") = R_NilValue,
            Rcpp::Named("nconv") = Rcpp::wrap(nconv),
            Rcpp::Named("niter") = Rcpp::wrap(niter)
        );
    }
    
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
    
    transformEigenvalues(evalsRm);
    
    return Rcpp::List::create(Rcpp::Named("values") = evalsRm,
                              Rcpp::Named("vectors") = eigenvectors,
                              Rcpp::Named("nconv") = Rcpp::wrap(nconv),
                              Rcpp::Named("niter") = Rcpp::wrap(niter));
}
