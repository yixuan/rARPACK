##' (Truncated) Singular Value Decomposition of a Matrix
##'
##' @description
##' Given an \code{m} by \code{n} matrix \code{A},
##' function \code{svds()} can find its largest \code{k}
##' singular values and the corresponding singular vectors.
##' It's also called the truncated singular value decomposition
##' since it only contains a subset of the whole singular triplets.
##' 
##' Currently \code{svds()} supports matrices of class "matrix",
##' "dgCMatrix" and "dsyMatrix". The latter two classes are defined
##' in the \pkg{Matrix} package, representing sparse matrix and
##' symmetric matrix respectively. Note that when \code{A} is symmetric,
##' SVD reduces to eigen decomposition, so you may consider using
##' \code{\link{eigs}()} instead.
##' 
##' @param A The matrix whose SVD is to be computed.
##' @param k Number of singular values requested.
##' @param nu Number of left singular vectors to be computed. This must
##'           be between 0 and \code{k}.
##' @param nv Number of right singular vectors to be computed. This must
##'           be between 0 and \code{k}.
##' @param opts Control parameters related to the computing
##'             algorithm. See Details below.
##' @param \dots Currently not used.
##'
##' @details The \code{opts} argument is a list that can supply any of the
##' following parameters:
##'
##' \describe{
##' \item{\code{ncv}}{Number of Lanzcos basis vectors to use. More vectors
##'                   will result in faster convergence, but with greater
##'                   memory use. \code{ncv} must be satisfy
##'                   \eqn{k+1 < ncv \le p}{k+1 < ncv <= p} where
##'                   \code{p = min(m, n)}.
##'                   Default is \code{min(p-1, max(2*k+1, 20))}.}
##' \item{\code{tol}}{Precision parameter. Default is 1e-10.}
##' \item{\code{maxitr}}{Maximum number of iterations. Default is 1000.}
##' }
##' 
##' @return A list with the following components:
##' \item{d}{A vector of the computed singular values.}
##' \item{u}{An \code{m} by \code{nu} matrix whose columns contain
##'          the left singular vectors. If \code{nu == 0}, \code{NULL}
##'          will be returned.}
##' \item{v}{An \code{n} by \code{nv} matrix whose columns contain
##'          the right singular vectors. If \code{nv == 0}, \code{NULL}
##'          will be returned.}
##' \item{nconv}{Number of converged singular values.}
##' \item{niter}{Number of iterations.}
##' @author Yixuan Qiu <\url{http://statr.me}>
##' @seealso \code{\link[base]{eigen}()}, \code{\link[base]{svd}()},
##' \code{\link[rARPACK]{eigs}()}.
##'
##' @export
##' @rdname svds
##' @keywords array
##' @examples
##' m = 100;
##' n = 20;
##' k = 5;
##' set.seed(111);
##' A = matrix(rnorm(m * n), m);
##' svds(A, k);
##' svds(A, k, nu = 0, nv = 3);
##'
##' ### more examples in examples/svds.R ###
svds <- function(A, k, nu = k, nv = k, opts = list(), ...)
    UseMethod("svds");

##' @rdname svds
##' @method svds matrix
##' @S3method svds matrix
##' @export
svds.matrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
    svds.real_gen(A, k, nu, nv, opts, ..., mattype = "matrix");

##' @rdname svds
##' @method svds dgeMatrix
##' @S3method svds dgeMatrix
##' @export
svds.dgeMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
    svds.real_gen(A, k, nu, nv, opts, ..., mattype = "dgeMatrix");

##' @rdname svds
##' @method svds dgCMatrix
##' @S3method svds dgCMatrix
##' @export
svds.dgCMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
    svds.real_gen(A, k, nu, nv, opts, ..., mattype = "dgCMatrix");

##' @rdname svds
##' @method svds dsyMatrix
##' @S3method svds dsyMatrix
##' @export
svds.dsyMatrix <- function(A, k, nu = k, nv = k, opts = list(), ...)
    svds.real_sym(A, k, nu, nv, opts, ..., mattype = "dsyMatrix");
