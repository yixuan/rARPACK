eigs.real_gen <- function(A, k, which, sigma, opts = list(), ...,
                          mattype = c("matrix", "dgeMatrix", "dgCMatrix", "dgRMatrix"))
{
    n = nrow(A);
    # Check whether 'A' is a square matrix
    if (n != ncol(A))
        stop("'A' must be a square matrix");
    # eigs() is not suitable for small matrices
    if (n < 3)
        stop("dimension of 'A' must be at least 3");
    
    # Matrix will be passed to C++, so we need to check the type.
    # ARPACK only supports matrices in float or double, so we need
    # to do the conversion if A is stored other than double.
    #
    # However, for sparse matrices defined in Matrix package,
    # they are always double, so we can omit this check. 
    if (mattype == "matrix" & typeof(A) != "double")
    {
        mode(A) = "double";
    }
    # Check the value of 'k'
    if (k <= 0 | k >= n - 1)
        stop("'k' must satisfy 0 < k < nrow(A) - 1.\nTo calculate all eigenvalues, try eigen()");
    
    # Check sigma
    # workmode == 1: ordinary
    # workmode == 3: Shift-invert mode
    if (is.null(sigma))
    {
        workmode = 1L;
        sigma = 0;
        sigmareal = TRUE;
    } else {
        workmode = 3L;
        if(abs(Im(sigma)) < 1e-17) sigma = Re(sigma);
        sigmareal = !is.complex(sigma);
    }
    
    # Arguments to be passed to ARPACK
    arpack.param = list(which = which,
                        ncv = min(n - 1, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000,
                        retvec = TRUE,
                        sigmar = Re(sigma[1]),
                        sigmai = Im(sigma[1]),
                        workmode = workmode);
    
    # Check the value of 'which'
    eigenv.type = c("LM", "SM", "LR", "SR", "LI", "SI");
    if (!(arpack.param$which %in% eigenv.type))
    {
        stop(sprintf("argument 'which' must be one of\n%s",
                     paste(eigenv.type, collapse = ", ")));
    }
    
    # Update parameters from 'opts' argument
    arpack.param[names(opts)] = opts;
    
    # Check the value of 'ncv'
    if (arpack.param$ncv < k + 2 | arpack.param$ncv > n)
        stop("'opts$ncv' must be >= k+2 and <= nrow(A)");
    
    # Call the C++ function
    res = .Call("eigs_gen",
                A,
                as.integer(n), as.integer(k),
                as.list(arpack.param),
                as.integer(MATTYPES[mattype]),
                PACKAGE = "rARPACK");
    
    # When workmode == 3 and sigmai != 0, we need to transform back
    # the eigenvalues
    # We only obtain
    #   nu = 0.5 * (1 / (lambda - sigma)) + 1 / (lambda - conj(sigma)))
    # So we need to solve lambda
    #   lambda = sigmar + (1 \pm sqrt(1 - 4 * sigmai^2 * nu^2)) / (2 * nu)
    # Use A * x = lambda * x to choose the correct root
    if(workmode == 3L & !sigmareal)
    {
        nu = res$values;
        sigmar = Re(sigma[1]);
        sigmai = Im(sigma[1]);
        tmp1 = sigmar + 0.5 / nu;
        tmp2 = 0.5 * sqrt(1 / nu^2 - 4 * sigmai^2);
        lambda1 = tmp1 + tmp2;
        lambda2 = tmp1 - tmp2;
        Ax.true = drop(as.numeric(A[1, ]) %*% res$vectors);
        lambdax1 = lambda1 * res$vectors[1, ];
        lambdax2 = lambda2 * res$vectors[1, ];
        mod1 = Mod(Ax.true - lambdax1);
        mod2 = Mod(Ax.true - lambdax2);
        res$values = ifelse(mod1 < mod2, lambda1, lambda2);
    }
    
    return(res);
}
