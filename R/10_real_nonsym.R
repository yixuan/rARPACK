eigs.real_nonsym <- function(A, k, which, sigma, opts = list(), ...,
                             mattype = c("dense", "sparse"))
{
    n = nrow(A);
    if (n != ncol(A))
        stop("'A' must be a square matrix");
    if (n < 3)
        stop("dimension of 'A' must be at least 3");
    
    # Matrix will be passed to C++, so we need to check the type.
    # ARPACK only supports matrices in float or double, so we need
    # to do the conversion if A is stored other than double.
    #
    # However, for sparse matrices defined in Matrix package,
    # they are always double, so we can omit this check. 
    if (mattype == "dense" & typeof(A) != "double")
    {
        mode(A) = "double";
    }
    if (k <= 0 | k >= n - 1)
        stop("'k' must satisfy 0 < k < nrow(A) - 1.\nTo calculate all eigenvalues, try eigen()");
    
    arpack.param = list(which = "LM",
                        ncv = min(n - 1, max(2 * k + 1, 20)),
                        tol = 1e-8,
                        maxitr = 300,
                        sigmar = 0.0,
                        sigmai = 0.0);
    eigenv.type = c("LM", "SM", "LR", "SR", "LI", "SI");
    if (inherits(sigma, "character"))
    {
        if (!(sigma %in% eigenv.type))
            stop(sprintf("when type of 'sigma' is character, it must one of\n%s",
                         paste(eigenv.type, collapse = ", ")));
        arpack.param$which = sigma;
    } else if (inherits(sigma, "numeric")){
        arpack.param$sigmar = Re(sigma[1]);
        arpack.param$sigmai = Im(sigma[1]);
    } else {
        stop("'sigma' must be either of the class 'character' or 'numeric'");
    }
    arpack.param[names(opts)] = opts;
    if (arpack.param$ncv < k + 2 | arpack.param$ncv > n)
        stop("'opts$ncv' must be >= k+2 and <= nrow(A)");
    
    funname = switch(mattype,
                     dense = "den_real_nonsym",
                     sparse = "sparse_real_nonsym",
                     stop("invalid value of 'mattype'"));
    
    res = .Call(funname, A, as.integer(n), as.integer(k),
                as.character(arpack.param$which), as.integer(arpack.param$ncv),
                as.numeric(arpack.param$tol), as.integer(arpack.param$maxitr),
                as.numeric(arpack.param$sigmar), as.numeric(arpack.param$sigmai),
                PACKAGE = "rarpack");
    
    if (!is.null(res$values))
    {
        res$values = res$values[1:res$nconv];
        return(res);
    } else {
        res$values = res$d.real[1:res$nconv] + res$d.imag[1:res$nconv] * 1i;
        res$vectors = res$v.real + res$v.imag * 1i;
        res$d.real = NULL;
        res$d.imag = NULL;
        res$v.real = NULL;
        res$v.imag = NULL;
    }
    return(res);
}
