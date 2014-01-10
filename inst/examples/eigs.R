library(rarpack);
n = 1000;
k = 5;

# test whether the calculated eigenvalues and eigenvectors satisfy
#                         A * x = lambda * x
eigen_resid = function(x, e)
{
    resid = x %*% e$vectors - e$vectors %*% diag(e$values);
    return(range(abs(resid)));
}

######################################
#
# Test for dense, non-symmetric matrix
#
######################################
set.seed(123);
x = matrix(rnorm(n^2), n);
x[sample(n^2, floor(n^2 / 2))] = 0;

# test whether the eigenvalues/eigenvectors are correct
# use default options
res1.1 = eigs(x, k);
eigen_resid(x, res1.1);
# which = "SM"
res1.2 = eigs(x, k, which = "SM");                            # ERROR
eigen_resid(x, res1.2);
# which = "LR"
res1.3 = eigs(x, k, which = "LR");
eigen_resid(x, res1.3);
# which = "SR"
res1.4 = eigs(x, k, which = "SR");
eigen_resid(x, res1.4);
# which = "LI"
res1.5 = eigs(x, k, which = "LI");
eigen_resid(x, res1.5);
# which = "SI"
res1.6 = eigs(x, k, which = "SI");
eigen_resid(x, res1.6);
# sigma = 0
res1.7 = eigs(x, k, sigma = 0);
eigen_resid(x, res1.7);
# sigma = 2
res1.8 = eigs(x, k, sigma = 2);
eigen_resid(x, res1.8);
# sigma = 1 + 1i
res1.9 = eigs(x, k, sigma = 1 + 1i);
eigen_resid(x, res1.9);

# only return eigenvalues
res1.10 = eigs(x, k, opts = list(retvec = FALSE));
range(abs(Mod(res1.10$values) - Mod(res1.1$values)));         # WRONG ORDER

######################################
#
# Test for sparse, non-symmetric matrix
#
######################################
xsp = as(x, "dgCMatrix");

# test whether the eigenvalues/eigenvectors are correct
# use default options
res2.1 = eigs(xsp, k);
eigen_resid(x, res2.1);
# which = "SM"
res2.2 = eigs(xsp, k, which = "SM");                            # ERROR
eigen_resid(x, res2.2);
# which = "LR"
res2.3 = eigs(xsp, k, which = "LR");
eigen_resid(x, res2.3);
# which = "SR"
res2.4 = eigs(xsp, k, which = "SR");
eigen_resid(x, res2.4);
# which = "LI"
res2.5 = eigs(xsp, k, which = "LI");
eigen_resid(x, res2.5);
# which = "SI"
res2.6 = eigs(xsp, k, which = "SI");
eigen_resid(x, res2.6);
# sigma = 0
res2.7 = eigs(xsp, k, sigma = 0);
eigen_resid(x, res2.7);
# sigma = 2
res2.8 = eigs(xsp, k, sigma = 2);
eigen_resid(x, res2.8);
# sigma = 1 + 1i
res2.9 = eigs(xsp, k, sigma = 1 + 1i);
eigen_resid(x, res2.9);

# only return eigenvalues
res2.10 = eigs(xsp, k, opts = list(retvec = FALSE));
range(abs(Mod(res2.10$values) - Mod(res2.1$values)));


######################################
#
# Test for dense, symmetric matrix
#
######################################
y = crossprod(x);
ysy = as(y, "dsyMatrix");

# test whether the eigenvalues/eigenvectors are correct
# use default options
res3.1 = eigs.sym(y, k);
res3.2 = eigs(ysy, k);
eigen_resid(y, res3.1);
eigen_resid(y, res3.2);
# which = "SM"
res3.3 = eigs.sym(y, k, which = "SM");                        # WARNING
res3.4 = eigs(ysy, k, which = "SM");                          # WARNING
eigen_resid(y, res3.3);
eigen_resid(y, res3.4);
# which = "LA"
res3.5 = eigs.sym(y, k, which = "LA");
res3.6 = eigs(ysy, k, which = "LA");
eigen_resid(y, res3.5);
eigen_resid(y, res3.6);
# which = "SA"
res3.7 = eigs.sym(y, k, which = "SA");                        # WARNING
res3.8 = eigs(ysy, k, which = "SA");                          # WARNING
eigen_resid(y, res3.7);
eigen_resid(y, res3.8);
# which = "BE"
res3.9 = eigs.sym(y, k, which = "BE");                        # WARNING
res3.10 = eigs(ysy, k, which = "BE");                         # WARNING
eigen_resid(y, res3.9);
eigen_resid(y, res3.10);
# sigma = 0
res3.11 = eigs.sym(y, k, sigma = 0);
res3.12 = eigs(ysy, k, sigma = 0);
eigen_resid(y, res3.11);
eigen_resid(y, res3.12);
# sigma = 2
res3.13 = eigs.sym(y, k, sigma = 2);
res3.14 = eigs(ysy, k, sigma = 2);
eigen_resid(y, res3.13);
eigen_resid(y, res3.14);

# only return eigenvalues
res3.15 = eigs.sym(y, k, opts = list(retvec = FALSE));
res3.16 = eigs(ysy, k, opts = list(retvec = FALSE));
range(abs(res3.15$values - res3.1$values));
range(abs(res3.16$values - res3.2$values));

