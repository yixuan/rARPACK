\name{eigs}
\alias{eigs}
\alias{eigs.sym}
\title{Find a few eigenvalues and eigenvectors for square matrix}
\usage{
  eigs(A, k, which = "LM", sigma = NULL, opts = list(),
    ...)

  eigs.sym(A, k, which = "LM", sigma = NULL, opts = list(),
    ..., lower = TRUE)
}
\arguments{
  \item{A}{The matrix whose eigen values/vectors are to be
  computed.}

  \item{k}{Number of eigenvalues requested.}

  \item{which}{Selection criteria. See Details below.}

  \item{sigma}{Shift parameter. See Details below.}

  \item{opts}{Control parameters related to the computing
  algorithm. See Details below.}

  \item{\dots}{Currently not used.}

  \item{lower}{For symmetric matrices, should the lower
  triangle or upper triangle be used.}
}
\value{
  A list of converged eigenvalues and eigenvectors.
  \item{nconv}{Number of converged eigenvalues.}
  \item{values}{Computed eigenvalues.}
  \item{vectors}{Computed eigenvectors. \code{vectors[, j]}
  corresponds to \code{values[j]}.}
}
\description{
  Given an \code{n} by \code{n} matrix \code{A}, function
  \code{eigs()} can calculate a limited number of
  eigenvalues and eigenvectors of \code{A}. Users can
  specify the selection criteria by argument \code{which},
  e.g., choosing the \code{k} largest or smallest
  eigenvalues and the corresponding eigenvectors.

  Currently \code{eigs()} supports matrices of class
  "matrix", "dgCMatrix" and "dsyMatrix".

  \code{eigs.sym()} assumes the matrix is symmetric, and
  only the lower triangle (or upper triangle, which is
  controlled by the argument \code{lower}) is used for
  computation, which in general reduces the workload.
  Notice that \code{eigs.sym()} only applies to "ordinary"
  matrix, i.e., of class "matrix". If you want to calculate
  eigen values/vectors of matrix of "dsyMatrix" class, use
  \code{eigs()} instead.
}
\details{
  The \code{which} argument is a character string that
  specifies the type of eigenvalues to be computed.
  Possible values are:

  \describe{ \item{"LM"}{The k eigenvalues with largest
  magnitude. Here the magnitude means the euclidean norm of
  complex numbers.} \item{"SM"}{The k eigenvalues with
  smallest magnitude. Here the magnitude means the
  euclidean norm of complex numbers.} \item{"LR"}{The k
  eigenvalues with largest real part.} \item{"SR"}{The k
  eigenvalues with smallest real part.} \item{"LI"}{The k
  eigenvalues with largest imaginary part.} \item{"SI"}{The
  k eigenvalues with smallest imaginary part.}
  \item{"LA"}{The k largest (algebraic) eigenvalues,
  considering any negative sign.} \item{"SA"}{The k
  smallest (algebraic) eigenvalues, considering any
  negative sign.} \item{"BE"}{Compute k eigenvalues, half
  from each end of the spectrum.  When k is odd, compute
  more from the high and then from the low end.} }

  \code{eigs()} with matrix type "matrix" and "dgCMatrix"
  can use "LM", "SM", "LR", "SR", "LI" and "SI".

  \code{eigs.sym()} and \code{eigs()} with matrix type
  "dsyMatrix" can use "LM", "SM", "LA", "SA" and "BE".

  The \code{sigma} argument is used in the Shift-and-invert
  mode. TODO

  The \code{opts} argument is a list that can supply any of
  the following parameters:

  \describe{ \item{\code{ncv}}{Number of Lanzcos basis
  vectors to use. More vectors will result in faster
  convergence, but with greater memory use. \code{ncv} must
  be satisfy \eqn{k+2\le ncv \le n}{k+2 <= ncv <= n}.
  Default is \code{min(n-1, max(2*k+1, 20))}.}
  \item{\code{tol}}{Precision parameter. Default is 1e-8.}
  \item{\code{maxitr}}{Maximum number of iterations.
  Default is 300.} \item{\code{retvec}}{Whether to compute
  eigenvectors. If FALSE, only calculate and return
  eigenvalues.} }
}
\examples{
n = 20;
k = 5;
## will have complex eigenvalues
A1 = matrix(rnorm(n^2), n);
eigs(A1, k);
## only have real eigenvalues,
## since A2 is symmetric
A2 = crossprod(A1);
eigs(A2, k);
eigs.sym(A2, k);
}
\author{
  Yixuan Qiu <\url{http://statr.me}>
}
\seealso{
  \code{\link[base]{eigen}()}, \code{\link[base]{svd}()}
}
\keyword{array}

