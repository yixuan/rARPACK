\docType{methods}
\name{eigs}
\alias{eigs}
\alias{eigs-matrix,numeric,character,missing,missing-method}
\alias{eigs-matrix,numeric,character,numeric,list-method}
\alias{eigs-matrix,numeric,character,numeric,missing-method}
\alias{eigs-matrix,numeric,missing,missing,missing-method}
\title{Find a few eigenvalues and eigenvectors}
\arguments{
  \item{A}{The matrix to compute eigenvalues from.}

  \item{k}{Number of eigenvalues requested.}

  \item{which}{Selection criteria. See Details below.}

  \item{sigma}{Shift parameter. See Details below.}

  \item{opts}{Control parameters related to the computing
  algorithm. See Details below.}

  \item{...}{Currently not used.}
}
\value{
  A list of converged eigenvalues and eigenvectors.
}
\description{
  Given an n by n matrix A, this function can calculate a
  limited number of eigenvalues and eigenvectors of A.
  Users can specify the selection criteria by argument
  \code{which}, e.g., choosing the k largest or smallest
  eigenvalues and the corresponding eigenvectors.
}
\examples{
A = crossprod(matrix(rnorm(10000), 100))
eigs(A, 6)
}
\seealso{
  \code{\link[base]{eigen}}, \code{\link[base]{svd}}
}
\keyword{array}
