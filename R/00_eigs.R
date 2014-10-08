##' Find a Specified Number of Eigenvalues/vectors for Square Matrix
##'
##' @description
##' Given an \code{n} by \code{n} matrix \code{A},
##' function \code{eigs()} can calculate a limited
##' number of eigenvalues and eigenvectors of \code{A}.
##' Users can specify the selection criteria by argument
##' \code{which}, e.g., choosing the \code{k} largest or smallest
##' eigenvalues and the corresponding eigenvectors.
##' 
##' Currently \code{eigs()} supports matrices of class "matrix"
##' (the most commonly used matrix type),
##' "dgeMatrix" (general matrix, equivalent to "matrix"),
##' "dgCMatrix" (sparse matrix), "dgRMatrix" (sparse matrix, row oriented)
##' and "dsyMatrix" (symmetric matrix).
##' All classes above except "matrix" are defined in the
##' \pkg{Matrix} package.
##' 
##' \code{eigs_sym()} assumes the matrix is symmetric,
##' and only the lower triangle (or upper triangle, which is
##' controlled by the argument \code{lower}) is used for
##' computation, which in some cases reduces the workload.
##' Notice that \code{eigs_sym()} only applies to "ordinary" matrix,
##' i.e., of class "matrix". If you want to calculate
##' eigen values/vectors of matrix of "dsyMatrix" class, use
##' \code{eigs()} instead.
##' 
##' @param A The matrix whose eigen values/vectors are to be computed.
##' @param k Number of eigenvalues requested.
##' @param which Selection criteria. See \sQuote{Details} below.
##' @param sigma Shift parameter. See \sQuote{Details} below.
##' @param opts Control parameters related to the computing
##' algorithm. See \sQuote{Details} below.
##' @param \dots Currently not used.
##' @param lower For symmetric matrices, should the lower triangle
##'              or upper triangle be used. 
##'
##' @details The \code{which} argument is a character string
##' that specifies the type of eigenvalues to be computed.
##' Possible values are:
##'
##' \describe{
##' \item{"LM"}{The k eigenvalues with largest magnitude. Here the
##'             magnitude means the euclidean norm of complex numbers.}
##' \item{"SM"}{The k eigenvalues with smallest magnitude. Here the
##'             magnitude means the euclidean norm of complex numbers.}
##' \item{"LR"}{The k eigenvalues with largest real part.}
##' \item{"SR"}{The k eigenvalues with smallest real part.}
##' \item{"LI"}{The k eigenvalues with largest imaginary part.}
##' \item{"SI"}{The k eigenvalues with smallest imaginary part.}
##' \item{"LA"}{The k largest (algebraic) eigenvalues, considering any
##'             negative sign.}
##' \item{"SA"}{The k smallest (algebraic) eigenvalues, considering any
##'             negative sign.}
##' \item{"BE"}{Compute k eigenvalues, half from each end of the spectrum.
##'             When k is odd, compute more from the high and then from the low end.}
##' }
##'
##' \code{eigs()} with matrix type "matrix", "dgeMatrix", "dgCMatrix"
##' and "dgRMatrix" can use "LM",
##' "SM", "LR", "SR", "LI" and "SI".
##' 
##' \code{eigs_sym()} and \code{eigs()} with matrix type "dsyMatrix"
##' can use "LM", "SM", "LA", "SA" and "BE".
##' 
##' The \code{sigma} argument is used in the shift-and-invert mode.
##' When \code{sigma} is not \code{NULL}, the selection criteria specified
##' by argument \code{which} will apply to
##' \eqn{1/(\lambda-\sigma)}{1/(lambda - sigma)}
##' where \eqn{\lambda}{lambda} are the eigenvalues of \eqn{A}{A}.
##' For example, if \eqn{A}{A} is positive definite and
##' \eqn{\sigma=0}{sigma = 0}, then \code{which = "LM"} will select the
##' largest values of \eqn{1/\lambda}{1/lambda}, which turns out to select
##' the smallest eigenvalues of \eqn{A}{A}. This method is preferable
##' to \code{which = "SM"} in that ARPACK is good at finding large
##' eigenvalues rather than finding small ones. More explanation of the
##' shift-and-invert mode can be found in the SciPy document,
##' \url{http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html}.
##' 
##' The \code{opts} argument is a list that can supply any of the
##' following parameters:
##'
##' \describe{
##' \item{\code{ncv}}{Number of Lanzcos basis vectors to use. More vectors
##'                   will result in faster convergence, but with greater
##'                   memory use. For general matrix, \code{ncv} must satisfy
##'                   \eqn{k+2\le ncv \le n}{k+2 <= ncv <= n}, and
##'                   for symmetric matrix, the constraint is
##'                   \eqn{k < ncv \le n}{k < ncv <= n}.
##'                   Default is \code{min(n, max(2*k+1, 20))}.}
##' \item{\code{tol}}{Precision parameter. Default is 1e-10.}
##' \item{\code{maxitr}}{Maximum number of iterations. Default is 1000.}
##' \item{\code{retvec}}{Whether to compute eigenvectors. If FALSE,
##'                      only calculate and return eigenvalues.}
##' }
##' 
##' @return A list of converged eigenvalues and eigenvectors.
##' \item{values}{Computed eigenvalues.}
##' \item{vectors}{Computed eigenvectors. \code{vectors[, j]} corresponds to \code{values[j]}.}
##' \item{nconv}{Number of converged eigenvalues.}
##' \item{niter}{Number of iterations in the computation.}
##' @author Yixuan Qiu <\url{http://statr.me}>
##' @seealso \code{\link[base]{eigen}()}, \code{\link[base]{svd}()},
##'          \code{\link[rARPACK]{svds}()}
##'
##' @export
##' @rdname eigs
##' @keywords array
##' @examples
##' n = 20;
##' k = 5;
##' 
##' ## Will have complex eigenvalues
##' set.seed(111);
##' A1 = matrix(rnorm(n^2), n);
##' eigs(A1, k);
##' 
##' ## Only have real eigenvalues,
##' ## since A2 is symmetric
##' A2 = crossprod(A1);
##' eigs(A2, k);
##' eigs_sym(A2, k);
##' 
##' ## Find the smallest (in absolute value) k eigenvalues of A2
##' eigs_sym(A2, k, which = "SM")
##' ## Another way to do this: use the sigma argument
##' eigs_sym(A2, k, sigma = 0)
##' ## The results should be the same,
##' ## but the latter method is far more stable on large matrices
##'
##' ### more examples in examples/eigs.R ###
eigs <- function(A, k, which = "LM", sigma = NULL, opts = list(), n = NULL, ...)
    UseMethod("eigs");

##' @rdname eigs
##' @export
eigs.matrix <- function(A, k, which = "LM", sigma = NULL,
                        opts = list(), ...)
    eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "matrix");

##' @rdname eigs
##' @export
eigs.dgeMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "dgeMatrix");

##' @rdname eigs
##' @export
eigs.dgCMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "dgCMatrix");

##' @rdname eigs
##' @export
eigs.dgRMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "dgRMatrix");

##' @rdname eigs
##' @export
eigs.dsyMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_sym(A, k, which, sigma, opts, ...,
                  mattype = "dsyMatrix", lower = (A@uplo == "L"))

##' @rdname eigs
##' @export
eigs.function <- function(FUN, k, which = "LM", sigma = NULL,
                           opts = list(), n, ...)
    eigs.fun(FUN, n, k, which, sigma, opts, ...,
                  mattype = "FUNCTION")



##' @rdname eigs
##' @usage eigs_sym(A, k, which = "LM", sigma = NULL, opts = list(),
##'   ..., lower = TRUE)
##' @export
eigs_sym <- function(A, k, which = "LM", sigma = NULL, opts = list(), ..., lower = TRUE)
{
    if(is.matrix(A))
    {
        eigs.real_sym(A, k, which, sigma, opts, ...,
                      mattype = "matrix", lower = lower);
    } else {
        stop("unsupported matrix type");
    }
}

# Matrix types
MATTYPES = c("matrix" = 0L, "dgeMatrix" = 1L, "dsyMatrix" = 2L,
             "dgCMatrix" = 3L, "dgRMatrix" = 4L, "FUNCTION" = 5L);
