library(rARPACK);
library(svd);
library(irlba);
library(microbenchmark);
library(plyr);
library(ggplot2);

n = 1000;
k = 20;
nu = 20;
nv = 20;

set.seed(123);
x = matrix(rnorm(n^2), n);
x[sample(n^2, floor(n^2 / 2))] = 0;
xsp = as(x, "dgCMatrix");
y = crossprod(x);
ysy = as(y, "dsyMatrix");

# svd package
f = function(v) as.numeric(xsp %*% v);
tf = function(v) as.numeric(crossprod(xsp, v));
extx = extmat(f, tf, nrow(xsp), ncol(xsp));

# irlba package
matmul = function(A, B, transpose = FALSE)
{
    if(transpose) as.numeric(crossprod(A, B)) else as.numeric(A %*% B);
}

######################################
#
# Compare eigs() with eigen()
#
######################################

### Eigenvalues only
op1.1 = microbenchmark(eigen = eigen(x, only.values = TRUE),
                       eigs = eigs(x, k, opts = list(retvec = FALSE)),
                       eigs.sp = eigs(xsp, k, opts = list(retvec = FALSE)));

op1.2 = microbenchmark(eigen = eigen(y),
                       eigs = eigs(y, k),
                       eigen.sym = eigen(y, symmetric = TRUE),
                       eigs_sym1 = eigs(ysy, k),
                       eigs_sym2 = eigs_sym(y, k));

######################################
#
# Compare svds() with svd()
#
######################################

### Singular values only ###
op2.1 = microbenchmark(svd0 = svd(x, nu = 0, nv = 0),
                       svds0 = svds(x, k, nu = 0, nv = 0),
                       svd = svd(x, nu, nv),
                       svds = svds(x, k, nu, nv),
                       propack.svd = propack.svd(x, k),
                       trlan.svd = trlan.svd(x, k),
                       irlba = irlba(x, nu, nv),
                       svds.sp0 = svds(xsp, k, nu = 0, nv = 0),
                       svds.sp = svds(xsp, k, nu, nv),
                       propack.svd.sp = propack.svd(extx, k),
                       trlan.svd.sp = trlan.svd(extx, k),
                       irlba.sp = irlba(xsp, nu, nv, matmul = matmul),
                       times = 20);
# write.csv(op2.1, "svd-Rblas-nosse.csv");
# write.csv(op2.1, "svd-OpenBLAS-nosse.csv");
# write.csv(op2.1, "svd-OpenBLAS-sse4.csv");

d1 = read.csv("svd-Rblas-nosse.csv");
d2 = read.csv("svd-OpenBLAS-nosse.csv");

dat = rbind(d1, d2);
dat$blas = rep(c("Default", "OpenBLAS"), c(nrow(d1), nrow(d2)));

dat.med = ddply(dat, c("expr", "blas"), summarize,
                medtime = median(time) / 1e6);
dat.med$expr = reorder(dat.med$expr, dat.med$medtime, max);
dat.med$blas = factor(dat.med$blas, levels = c("OpenBLAS", "Default"))
ggplot(dat.med, aes(x = expr, y = medtime, fill = blas)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_hue("BLAS", limits = c("Default", "OpenBLAS")) +
    scale_y_continuous("Elapsed time (milliseconds)") +
    scale_x_discrete("Functions", labels = c("svds[sparse, value-only]",
                                             "svds[sparse]",
                                             "svds[dense, value-only]",
                                             "svds[dense]",
                                             "propack.svd[dense]",
                                             "trlan.svd[dense]",
                                             "svd[dense, value-only]",
                                             "propack.svd[sparse]",
                                             "irlba[dense]",
                                             "irlba[sparse]",
                                             "trlan.svd[sparse]",
                                             "svd[dense]")) +
    coord_flip() +
    theme_grey(base_size = 18);
