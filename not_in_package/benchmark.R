library(Matrix)
library(rARPACK)
library(svd)
library(irlba)
library(microbenchmark)
library(dplyr)
library(ggplot2)

n = 1000
k = 20
nu = 20
nv = 20

set.seed(123)
x = matrix(rnorm(n^2), n)
x[sample(n^2, floor(n^2 / 2))] = 0
xsp = as(x, "dgCMatrix")
y = crossprod(x)
ysy = as(y, "dsyMatrix")

## svd package
f = function(v) as.numeric(xsp %*% v)
tf = function(v) as.numeric(crossprod(xsp, v))
extx = extmat(f, tf, nrow(xsp), ncol(xsp))

## irlba package
matmul = function(A, B, transpose = FALSE)
{
    if(transpose) as.numeric(crossprod(A, B)) else as.numeric(A %*% B);
}



## Visualize results
visualize = function(res)
{
    dat = res %>% group_by(expr) %>%
        summarize(medtime = median(time) / 1e6) %>%
        mutate(expr = reorder(expr, medtime))
    
    g = ggplot(dat, aes(x = expr, y = medtime)) +
        geom_bar(stat = "identity") +
        scale_y_continuous("Elapsed time (milliseconds)") +
        scale_x_discrete("Functions") +
        coord_flip() +
        theme_grey(base_size = 18)
    print(g)
}

## How many times to run each function
nrun = 10

#########################################
#
# Compare eigs() with eigen()
#
#########################################

### Non-symmetric matrix
res1.1 = microbenchmark(
    "eigen[value-only]"        = eigen(x, only.values = TRUE),
    "eigs[value-only]"         = eigs(x, k, opts = list(retvec = FALSE)),
    "eigs[sparse, value-only]" = eigs(xsp, k, opts = list(retvec = FALSE)),
    "eigen"                    = eigen(x),
    "eigs"                     = eigs(x, k),
    "eigs[sparse]"             = eigs(xsp, k),
    times = nrun
)

visualize(res1.1)

### Symmetric matrix
res1.2 = microbenchmark(
    "eigen[value-only]"            = eigen(y, symmetric = TRUE, only.values = TRUE),
    "eigs[interface1, value-only]" = eigs_sym(y, k, opts = list(retvec = FALSE)),
    "eigs[interface2, value-only]" = eigs(ysy, k, opts = list(retvec = FALSE)),
    "eigen"                        = eigen(y, symmetric = TRUE),
    "eigs[interface1]"             = eigs_sym(y, k),
    "eigs[interface2]"             = eigs(ysy, k),
    times = nrun
)

visualize(res1.2)

#########################################
#
# Compare svds() with other SVD functions
#
#########################################

res2.1 = microbenchmark(
    "svd[value-only]"          = svd(x, nu = 0, nv = 0),
    "svds[value-only]"         = svds(x, k, nu = 0, nv = 0),
    "svd"                      = svd(x, nu, nv),
    "svds"                     = svds(x, k, nu, nv),
    "propack"                  = propack.svd(x, k),
    "trlan"                    = trlan.svd(x, k),
    "irlba"                    = irlba(x, nu, nv),
    "svds[sparse, value-only]" = svds(xsp, k, nu = 0, nv = 0),
    "svds[sparse]"             = svds(xsp, k, nu, nv),
    "propack[sparse]"          = propack.svd(extx, k),
    "trlan[sparse]"            = trlan.svd(extx, k),
    "irlba[sparse]"            = irlba(xsp, nu, nv, matmul = matmul),
    times = nrun
)

visualize(res2.1)
