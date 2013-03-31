##' @rdname eigs-methods
##' @aliases eigs,dgCMatrix,numeric,character,numeric,list-method
setMethod("eigs", signature(A = "dgCMatrix",
                            k = "numeric",
                            which = "character",
                            sigma = "numeric",
                            opts = "list"),
          function(A, k, which, sigma, opts, ...)
              eigs.real_nonsym(A, k, which, sigma, opts, ...,
                               mattype = "sparse"));
##' @rdname eigs-methods
##' @aliases eigs,dgCMatrix,numeric,character,numeric,missing-method
setMethod("eigs", signature(A = "dgCMatrix",
                            k = "numeric",
                            which = "character",
                            sigma = "numeric",
                            opts = "missing"),
          function(A, k, which, sigma, ...)
              eigs.real_nonsym(A, k, which, sigma, list(), ...,
                               mattype = "sparse"));
##' @rdname eigs-methods
##' @aliases eigs,dgCMatrix,numeric,character,missing,missing-method
setMethod("eigs", signature(A = "dgCMatrix",
                            k = "numeric",
                            which = "character",
                            sigma = "missing",
                            opts = "missing"),
          function(A, k, which, ...)
              eigs.real_nonsym(A, k, which, 0.0, list(), ...,
                               mattype = "sparse"));
##' @rdname eigs-methods
##' @aliases eigs,dgCMatrix,numeric,missing,missing,missing-method
setMethod("eigs", signature(A = "dgCMatrix",
                            k = "numeric",
                            which = "missing",
                            sigma = "missing",
                            opts = "missing"),
          function(A, k, which, ...)
              eigs.real_nonsym(A, k, "LM", 0.0, list(), ...,
                               mattype = "sparse"));
