svds.real_nonsym <- function(A, k, nu = k, nv = k, opts = list(), ...,
                             mattype = c("matrix", "dgCMatrix"))
{
    m = nrow(A);
    n = ncol(A);
    wd = min(m, n);
    
    # Check for matrices that are too small
    if (wd < 3)
        stop("nrow(A) and ncol(A) should be at least 3");
    
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
    if (k <= 0 | k >= wd - 1)
        stop("'k' must satisfy 0 < k < min(nrow(A), ncol(A)) - 1.\nTo calculate all singular values, try svd()");
    
    # Check the values of 'nu' and 'nv'
    if (nu < 0 | nv < 0 | nu > k | nv > k)
        stop("'nu' and 'nv' must satisfy 0 <= nu <= k and 0 <= nv <= k");
    
    # Arguments to be passed to ARPACK
    arpack.param = list(ncv = min(wd - 1, max(2 * k + 1, 20)),
                        tol = 1e-10,
                        maxitr = 1000);
    
    # Update parameters from 'opts' argument
    arpack.param[names(opts)] = opts;
    
    # Check the value of 'ncv'
    if (arpack.param$ncv < k + 2 | arpack.param$ncv > wd)
        stop("'opts$ncv' must be >= k+2 and <= min(nrow(A), ncol(A))");
    
    # Different names of calls according to the type of matrix
    funname = switch(mattype,
                     matrix = "den_real_nonsym_svd",
                     dgCMatrix = "sparse_real_nonsym_svd",
                     stop("invalid matrix type"));
    
    # Calling the C++ function
    res = .Call(funname,
                A,
                as.integer(m), as.integer(n),
                as.integer(k), as.integer(nu), as.integer(nv),
                as.list(arpack.param),
                PACKAGE = "rARPACK");
    
    return(res);
}
