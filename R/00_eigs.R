##' Find a Specified Number of Eigenvalues/vectors for Square Matrix
##'
##' @description
##' Given an \eqn{n} by \eqn{n} matrix \eqn{A},
##' function \code{eigs()} can calculate a limited
##' number of eigenvalues and eigenvectors of \eqn{A}.
##' Users can specify the selection criteria by argument
##' \code{which}, e.g., choosing the \eqn{k} largest or smallest
##' eigenvalues and the corresponding eigenvectors.
##' 
##' Currently \code{eigs()} supports matrices of the following classes:
##' 
##' \tabular{ll}{
##'   \code{matrix}     \tab The most commonly used matrix type,
##'                          defined in \strong{base} package.\cr
##'   \code{dgeMatrix}  \tab General matrix, equivalent to \code{matrix},
##'                          defined in \strong{Matrix} package.\cr
##'   \code{dgCMatrix}  \tab Column oriented sparse matrix, defined in
##'                          \strong{Matrix} package.\cr
##'   \code{dgRMatrix}  \tab Row oriented sparse matrix, defined in
##'                          \strong{Matrix} package.\cr
##'   \code{dsyMatrix}  \tab Symmetrix matrix, defined in \strong{Matrix}
##'                          package.
##' }
##' 
##' \code{eigs_sym()} assumes the matrix is symmetric,
##' and only the lower triangle (or upper triangle, which is
##' controlled by the argument \code{lower}) is used for
##' computation, which in some cases reduces the workload.
##' Notice that \code{eigs_sym()} only applies to "ordinary" matrix,
##' i.e., of class "matrix". If you want to calculate
##' eigenvalues/vectors of matrix of "dsyMatrix" class, use
##' \code{eigs()} instead.
##' 
##' @param A The matrix whose eigenvalues/vectors are to be computed.
##' @param k Number of eigenvalues requested.
##' @param which Selection criteria. See \strong{Details} below.
##' @param sigma Shift parameter. See section \strong{Shift-And-Invert Mode}.
##' @param opts Control parameters related to the computing
##'             algorithm. See \strong{Details} below.
##' @param \dots Currently not used.
##' @param lower For symmetric matrices, should the lower triangle
##'              or upper triangle be used.
##' @param args Argument passed to \code{A} when \code{A} is a
##'             function. See section \strong{Function Interface}
##'             for details.
##'
##' @details The \code{which} argument is a character string
##' that specifies the type of eigenvalues to be computed.
##' Possible values are:
##'
##' \tabular{ll}{
##'   "LM"  \tab  The \eqn{k} eigenvalues with largest magnitude. Here the
##'               magnitude means the Euclidean norm of complex numbers.\cr
##'   "SM"  \tab  The \eqn{k} eigenvalues with smallest magnitude.\cr
##'   "LR"  \tab  The \eqn{k} eigenvalues with largest real part.\cr
##'   "SR"  \tab  The \eqn{k} eigenvalues with smallest real part.\cr
##'   "LI"  \tab  The \eqn{k} eigenvalues with largest imaginary part.\cr
##'   "SI"  \tab  The \eqn{k} eigenvalues with smallest imaginary part.\cr
##'   "LA"  \tab  The \eqn{k} largest (algebraic) eigenvalues, considering any
##'               negative sign.\cr
##'   "SA"  \tab  The \eqn{k} smallest (algebraic) eigenvalues, considering any
##'               negative sign.\cr
##'   "BE"  \tab  Compute \eqn{k} eigenvalues, half from each end of the
##'               spectrum. When \eqn{k} is odd, compute more from the high
##'               and then from the low end.
##' }
##'
##' \code{eigs()} with matrix type "matrix", "dgeMatrix", "dgCMatrix"
##' and "dgRMatrix" can use "LM",
##' "SM", "LR", "SR", "LI" and "SI".
##' 
##' \code{eigs_sym()}, and \code{eigs()} with matrix type "dsyMatrix"
##' can use "LM", "SM", "LA", "SA" and "BE".
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
##' @section Shift-And-Invert Mode:
##' The \code{sigma} argument is used in the shift-and-invert mode.
##' 
##' When \code{sigma} is not \code{NULL}, the selection criteria specified
##' by argument \code{which} will apply to
##' 
##' \deqn{\frac{1}{\lambda-\sigma}}{1/(\lambda-\sigma)}
##' 
##' where \eqn{\lambda}'s are the eigenvalues of \eqn{A}. This mode is useful
##' when user wants to find eigenvalues closest to a given number.
##' For example, if \eqn{\sigma=0}, then \code{which = "LM"} will select the
##' largest values of \eqn{1/|\lambda|}, which turns out to select
##' eigenvalues of \eqn{A} that have the smallest magnitude. The result of
##' using \code{which = "LM", sigma = 0} will be the same as
##' \code{which = "SM"}, but the former one is preferable
##' in that ARPACK is good at finding large
##' eigenvalues rather than small ones. More explanation of the
##' shift-and-invert mode can be found in the SciPy document,
##' \url{http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html}.
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
##' n = 20
##' k = 5
##' 
##' ## Will have complex eigenvalues
##' set.seed(111)
##' A1 = matrix(rnorm(n^2), n)
##' eigs(A1, k)
##' 
##' f = function(x, args)
##' {
##'     mat = args$mat
##'     mat %*% x
##' }
##' eigs(f, k, args = list(n = n, mat = A1))
##' 
##' ## Only have real eigenvalues,
##' ## since A2 is symmetric
##' A2 = crossprod(A1)
##' eigs(A2, k)
##' eigs_sym(A2, k)
##' 
##' ## Find the smallest (in absolute value) k eigenvalues of A2
##' eigs_sym(A2, k, which = "SM")
##' ## Another way to do this: use the sigma argument
##' eigs_sym(A2, k, sigma = 0)
##' ## The results should be the same,
##' ## but the latter method is far more stable on large matrices
##'
##' ### more examples in examples/eigs.R ###
eigs <- function(A, k, which = "LM", sigma = NULL,
                 opts = list(), args = list(n = NULL), ...)
    UseMethod("eigs")

##' @rdname eigs
##' @export
eigs.matrix <- function(A, k, which = "LM", sigma = NULL,
                        opts = list(), ...)
{
    if(isSymmetric(A) &
           which %in% c("LM", "SM", "LR", "SR") &
           (is.null(sigma) || Im(sigma) == 0))
    {
        if(which == "LR")  which = "LA"
        if(which == "SR")  which = "SA"
        eigs.real_sym(A, k, which, sigma, opts, ..., mattype = "symmatrix")
    } else {
        eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "matrix")
    }
}

##' @rdname eigs
##' @export
eigs.dgeMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
{
    if(isSymmetric(A) &
           which %in% c("LM", "SM", "LR", "SR") &
           (is.null(sigma) || Im(sigma) == 0))
    {
        if(which == "LR")  which = "LA"
        if(which == "SR")  which = "SA"
        eigs.real_sym(A, k, which, sigma, opts, ..., mattype = "dgeMatrix")
    } else {
        eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "dgeMatrix")
    }
}

##' @rdname eigs
##' @export
eigs.dgCMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
{
    if(isSymmetric(A) &
           which %in% c("LM", "SM", "LR", "SR") &
           (is.null(sigma) || Im(sigma) == 0))
    {
        if(which == "LR")  which = "LA"
        if(which == "SR")  which = "SA"
        eigs.real_sym(A, k, which, sigma, opts, ..., mattype = "dgCMatrix")
    } else {
        eigs.real_gen(A, k, which, sigma, opts, ..., mattype = "dgCMatrix")
    }
}

##' @rdname eigs
##' @export
## isSymmetric() does not support dgRMatrix
eigs.dgRMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_gen(A, k, which, sigma, opts, ...,
                  mattype = "dgRMatrix")

##' @rdname eigs
##' @export
eigs.dsyMatrix <- function(A, k, which = "LM", sigma = NULL,
                           opts = list(), ...)
    eigs.real_sym(A, k, which, sigma, opts, ...,
                  mattype = "dsyMatrix", lower = (A@uplo == "L"))

##' @rdname eigs
##' @export
eigs.function <- function(A, k, which = "LM", sigma = NULL,
                          opts = list(), ..., args = list(n = NULL))
    eigs.fun(A, k, which, sigma, opts, ..., mattype = "function",
             args = args)



##' @rdname eigs
##' @usage eigs_sym(A, k, which = "LM", sigma = NULL, opts = list(),
##'   ..., lower = TRUE)
##' @export
eigs_sym <- function(A, k, which = "LM", sigma = NULL, opts = list(), ...,
                     lower = TRUE)
{
    if(is.matrix(A))
    {
        eigs.real_sym(A, k, which, sigma, opts, ...,
                      mattype = "matrix", lower = lower)
    } else {
        stop("unsupported matrix type")
    }
}

# Matrix types
MATTYPES = c("matrix" = 0L, "symmatrix" = 1L, "dgeMatrix" = 2L,
             "dsyMatrix" = 3L, "dgCMatrix" = 4L, "dgRMatrix" = 5L,
             "function" = 6L)
