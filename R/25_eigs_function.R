eigs.fun <- function(FUN, k, which, sigma, opts = list(),
                     args = list(n = NULL), ..., mattype = c("function"))
{
    # Check whether n is NULL
    n = args[["n"]]
    if (is.null(n))
        stop("must provide 'n' in 'args', the dimension of the implicit matrix")
    # Check whether FUN(x) returns n elements
    if (length(FUN(rnorm(n), args)) != n)
        stop("the provided function must return a length-n vector")
    # eigs() is not suitable for small matrices
    if (n < 3)
        stop("'n' must be at least 3")
    
    # Check the value of 'k'
    if (k <= 0 | k >= n)
        stop("'k' must satisfy 0 < k < n")
    
    # Check sigma, we do not have shift-and-invert mode for functional A
    if (!is.null(sigma))
    {
        warning("'sigma' is ignored when 'A' is a function")
    }
    
    # Arguments to be passed to ARPACK
    arpack.param = list(which = which,
                        ncv = min(n, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000,
                        retvec = TRUE,
                        sigma = 0,
                        workmode = 1)
    
    # Check the value of 'which'
    eigenv.type = c("LM", "SM", "LR", "SR", "LI", "SI")
    if (!(arpack.param$which %in% eigenv.type))
    {
        stop(sprintf("argument 'which' must be one of\n%s",
                     paste(eigenv.type, collapse = ", ")))
    }
    
    # Update parameters from 'opts' argument
    arpack.param[names(opts)] = opts
    
    # Check the value of 'ncv'
    if (arpack.param$ncv <= k | arpack.param$ncv > n)
        stop("'opts$ncv' must be > k and <= n")
    
    # Call the C++ function
    res = .Call("eigs_fun",
                as.function(FUN), as.list(args),
                as.integer(n), as.integer(k),
                as.list(arpack.param),
                as.integer(MATTYPES[mattype]),
                PACKAGE = "rARPACK")
    
    return(res)
}
