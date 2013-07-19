##' Find a few eigenvalues and eigenvectors for square matrix
##'
##' Given an \code{n} by \code{n} matrix \code{A},
##' this function can calculate a limited
##' number of eigenvalues and eigenvectors of \code{A}. Users can specify
##' the selection criteria by argument \code{which},
##' e.g., choosing the \code{k} largest or smallest
##' eigenvalues and the corresponding eigenvectors.
##' 
##' @param A The matrix whose eigen values/vectors are to be computed.
##' @param k Number of eigenvalues requested.
##' @param which Selection criteria. See Details below.
##' @param sigma Shift parameter. See Details below.
##' @param opts Control parameters related to the computing
##' algorithm. See Details below.
##' @param \dots Currently not used.
##'
##' @details The \code{which} argument is a character string
##' that specifies the type of eigenvalues to be computed.
##' Possible values are:
##'
##' \describe{
##' \item{"LM"}{The k eigenvalues with largest magnitude.}
##' \item{"SM"}{The k eigenvalues with smallest magnitude.}
##' \item{"LR"}{The k eigenvalues with largest real part.}
##' \item{"SR"}{The k eigenvalues with smallest real part.}
##' \item{"LI"}{The k eigenvalues with largest imaginary part.}
##' \item{"SI"}{The k eigenvalues with smallest imaginary part.}
##' }
##'
##' The \code{opts} argument is a list that can supply any of the
##' following parameters:
##'
##' \describe{
##' \item{\code{ncv}}{Number of Lanzcos basis vectors to use. More vectors will result in faster convergence, but with greater memory use. \code{ncv} must be satisfy \eqn{k+2\le ncv \le n}{k+2 <= ncv <= n}. Default is \code{min(n-1, max(2*k+1, 20))}.}
##' \item{\code{tol}}{Precision parameter. Default is 1e-8.}
##' \item{\code{maxitr}}{Maximum number of iterations. Default is 300.}
##' \item{\code{retvec}}{Whether to compute eigenvectors. If FALSE,
##'                      only returning eigenvalues.}
##' }
##' 
##' @return A list of converged eigenvalues and eigenvectors.
##' \item{nconv}{Number of converged eigenvalues.}
##' \item{values}{Computed eigenvalues.}
##' \item{vectors}{Computed eigenvectors. \code{vectors[, j]} corresponds to \code{values[j]}.}
##' @note Currently only real nonsymmetric matrices are supported.
##' @author Yixuan Qiu <\url{http://statr.me}>
##' @seealso \code{\link[base]{eigen}()}, \code{\link[base]{svd}()}
##'
##' @export
##' @rdname eigs
##' @keywords array
##' @examples
##' n = 20;
##' k = 5;
##' ## will have complex eigenvalues
##' A1 = matrix(rnorm(n^2), n);
##' eigs(A1, k);
##' ## only have real eigenvalues,
##' ## since A2 is symmetric
##' A2 = crossprod(A1);
##' eigs(A2, k);
eigs <- function(A, k, which = "LM", sigma = 0.0, opts = list(), ...)
{
    if(inherits(A, "matrix"))
    {
        eigs.real_nonsym(A, k, which, sigma, opts, ...,
                         mattype = "matrix");
    } else if(inherits(A, "dgCMatrix")) {
        eigs.real_nonsym(A, k, which, sigma, opts, ...,
                         mattype = "dgCMatrix");
    } else if(inherits(A, "dsyMatrix")){
        eigs.real_sym(A, k, which, sigma, opts, ...,
                      mattype = "dsyMatrix")
    } else {
        stop("unsupported matrix type");
    }
}

##' @rdname eigs
##' @export
eigs.sym <- function(A, k, which = "LM", sigma = 0.0, opts = list(), ..., lower = TRUE)
{
    if(inherits(A, "matrix"))
    {
        eigs.real_sym(A, k, which, sigma, opts, ...,
                      mattype = "matrix", lower = lower);
    } else {
        stop("unsupported matrix type");
    }
}
