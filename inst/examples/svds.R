library(rARPACK);
n = 100;
p = 50;
k = 5;

# test whether the calculated singular values and singular vectors
# are consistent with those calculated by svd()
svd_resid = function(rsvd, rsvds)
{
    k = length(rsvds$d);
    nu = ncol(rsvds$u);
    nv = ncol(rsvds$v);
    
    dres = abs(rsvd$d[1:k] - rsvds$d);
    ures = if(!is.null(nu))
               abs(rsvd$u[, 1:nu] %*% diag(sign(rsvd$u[1, 1:nu])) -
                   rsvds$u %*% diag(sign(rsvds$u[1, ]))) else NULL;
    vres = if(!is.null(nv))
               abs(rsvd$v[, 1:nv] %*% diag(sign(rsvd$v[1, 1:nv])) -
                   rsvds$v %*% diag(sign(rsvds$v[1, ]))) else NULL;
    cat("d:", range(dres), "\n");
    cat("u:", if(is.null(ures)) "unknown" else range(ures), "\n");
    cat("v:", if(is.null(vres)) "unknown" else range(vres), "\n");
}

######################################
#
# Test for dense, non-symmetric matrix
#
######################################
set.seed(123);
x = matrix(rnorm(n * p), n);
x[sample(n * p, floor(n * p / 2))] = 0;
xt = t(x);

# results by svd()
res1.0 = svd(x);
res1.0t = svd(xt);

# test whether the results are consistent with svd()
# use default options
res1.1 = svds(x, k);
res1.2 = svds(xt, k);
svd_resid(res1.0, res1.1);
svd_resid(res1.0t, res1.2);

# do not return U
res1.3 = svds(x, k, nu = 0);
res1.4 = svds(xt, k, nu = 0);
svd_resid(res1.0, res1.3);
svd_resid(res1.0t, res1.4);

# do not return V
res1.5 = svds(x, k, nv = 0);
res1.6 = svds(xt, k, nv = 0);
svd_resid(res1.0, res1.5);
svd_resid(res1.0t, res1.6);

# only singular values
res1.7 = svds(x, k, nu = 0, nv = 0);
res1.8 = svds(xt, k, nu  = 0, nv = 0);
svd_resid(res1.0, res1.7);
svd_resid(res1.0t, res1.8);

######################################
#
# Test for sparse, non-symmetric matrix
#
######################################
xsp = as(x, "dgCMatrix");
xspt = t(xsp);

# test whether the results are consistent with svd()
# use default options
res2.1 = svds(xsp, k);
res2.2 = svds(xspt, k);
svd_resid(res1.0, res2.1);
svd_resid(res1.0t, res2.2);

# do not return U
res2.3 = svds(xsp, k, nu = 0);
res2.4 = svds(xspt, k, nu = 0);
svd_resid(res1.0, res2.3);
svd_resid(res1.0t, res2.4);

# do not return V
res2.5 = svds(xsp, k, nv = 0);
res2.6 = svds(xspt, k, nv = 0);
svd_resid(res1.0, res2.5);
svd_resid(res1.0t, res2.6);

# only singular values
res2.7 = svds(xsp, k, nu = 0, nv = 0);
res2.8 = svds(xspt, k, nu  = 0, nv = 0);
svd_resid(res1.0, res2.7);
svd_resid(res1.0t, res2.8);


######################################
#
# Test for dense, symmetric matrix
#
######################################
y = crossprod(x);
ysy = as(y, "dsyMatrix");

# results by svd()
res3.0 = svd(y);

# test whether the results are consistent with svd()
# use default options
res3.1 = svds(ysy, k);
svd_resid(res3.0, res3.1);

# do not return U
res3.2 = svds(ysy, k, nu = 0);
svd_resid(res3.0, res3.2);

# do not return V
res3.3 = svds(ysy, k, nv = 0);
svd_resid(res3.0, res3.3);

# only singular values
res3.4 = svds(ysy, k, nu = 0, nv = 0);
svd_resid(res3.0, res3.4);

