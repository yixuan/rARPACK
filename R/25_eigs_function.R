eigs.fun <- function(FUN, n, k, which, sigma, opts = list(), ...,
                          mattype = c("function"))
{
    # Check whether n is NULL
    if (is.null(n))
        stop("Must provide the dimension of the implicit matrix")
    # Check whether 'A' is a square matrix
    if (length(FUN(rnorm(n))) != n)
        stop("The implicit matrix must return a length-n vector");
    # eigs() is not suitable for small matrices
    if (n < 3)
        stop("dimension of 'A' must be at least 3");
    
    # Check the value of 'k'
    if (k <= 0 | k >= n)
        stop("'k' must satisfy 0 < k < nrow(A)");
    
    # Check sigma
    # workmode == 1: ordinary
    # workmode == 3: Shift-invert mode
    if (is.null(sigma))
    {
        workmode = 1L;
        sigma = 0;
    } else {
        workmode = 3L;
        if(is.complex(sigma)) warning("only real part of sigma is used");
        sigma = Re(sigma);
    }
    
    # Arguments to be passed to ARPACK
    arpack.param = list(which = which,
                        ncv = min(n, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000,
                        retvec = TRUE,
                        sigma = sigma,
                        workmode = workmode);
    
    # Check the value of 'which'
    eigenv.type = c("LM", "SM", "LA", "SA", "BE");
    if (!(arpack.param$which %in% eigenv.type))
    {
        stop(sprintf("argument 'which' must be one of\n%s",
                     paste(eigenv.type, collapse = ", ")));
    }
    
    # Update parameters from 'opts' argument
    arpack.param[names(opts)] = opts;
    
    # Check the value of 'ncv'
    if (arpack.param$ncv <= k | arpack.param$ncv > n)
        stop("'opts$ncv' must be > k and <= nrow(A)");
    
    # Call the C++ function
    res = .Call("eigs_fun",
                FUN,
                as.integer(n), as.integer(k),
                as.list(arpack.param),
                as.integer(MATTYPES[mattype]),
                PACKAGE = "rARPACK");
    
    return(res);
}
