##' Find a few eigenvalues and eigenvectors
##'
##' Given an n by n matrix A, this function can calculate a limited
##' number of eigenvalues and eigenvectors of A. Users can specify
##' the selection criteria by argument \code{which},
##' e.g., choosing the k largest or smallest
##' eigenvalues and the corresponding eigenvectors.
##'
##' @param A The matrix to compute eigenvalues from.
##' @param k Number of eigenvalues requested.
##' @param which Selection criteria. See Details below.
##' @param sigma Shift parameter. See Details below.
##' @param opts Control parameters related to the computing
##' algorithm. See Details below.
##' @param ... Currently not used.
##'
##' @export
##' @docType methods
##' @rdname eigs-methods
##' 
##' @usage eigs(A, k, which, sigma, opts, ...)
##'
##' @return A list of converged eigenvalues and eigenvectors.
##'
##' @seealso \code{\link[base]{eigen}}, \code{\link[base]{svd}}
##' @keywords array
##' @examples
##' A = crossprod(matrix(rnorm(10000), 100));
##' eigs(A, 6);
setGeneric("eigs",
    function(A, k, which, sigma, opts, ...) standardGeneric("eigs"));
