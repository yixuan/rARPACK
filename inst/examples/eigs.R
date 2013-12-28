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










# test dense non-symmetric matrix
x = matrix(rnorm(n^2), n);
x[sample(n^2, floor(n^2 / 2))] = 0;
system.time(res1 <- eigs(x, k));
res1$values;
system.time(res1.5 <- eigs(x, k, sigma = 1));
res1.5$values;

# test sparse non-symmetric matrix
xsp = as(x, "dgCMatrix");
system.time(res2 <- eigs(xsp, k));
res2$values;
system.time(res2.5 <- eigs(xsp, k, sigma = 1));
res2.5$values;

# test dense symmetric matrix of class "matrix"
y = crossprod(x);
system.time(res3 <- eigs.sym(y, k));
system.time(res4 <- eigs.sym(y, k, lower = FALSE));
res3$values;
res4$values;
y[1, 5] = 2;
system.time(res5 <- eigs.sym(y, k));
system.time(res6 <- eigs.sym(y, k, lower = FALSE));
res5$values;
res6$values;

# test dense symmetric matrix of calss "dsyMatrix"
y = crossprod(x);
ysy = as(y, "dsyMatrix");
system.time(res7 <- eigs(ysy, k));
ysy2 = t(ysy);
system.time(res8 <- eigs(ysy2, k));
res7$values;
res8$values;



# test selection criteria
vec = eigen(crossprod(matrix(rnorm(36), 6, 6)))$vectors;
y = vec %*% diag(c(-3, -2, -1, 1, 2, 3)) %*% t(vec);
eigs(y, 3, "LM");
eigs(y, 3, "SM");
eigs(y, 3, "LR");
eigs(y, 3, "SR");
eigs(y, 3, "LI");
eigs(y, 3, "SI");
eigs(y, 3, "LM", sigma = 1.9);
eigs(y, 3, "LM", sigma = 2);
eigs(y, 3, "LM", sigma = 0);
eigs(y, 3, "LM", sigma = 1.9 + 0i);
eigs(y, 3, "LM", sigma = 1.9 + 1i); # Caution!
eigs.sym(y, 3, "LM");
eigs.sym(y, 3, "SM");
eigs.sym(y, 3, "LA");
eigs.sym(y, 3, "SA");
eigs.sym(y, 3, "BE");
eigs.sym(y, 3, "LM", sigma = 1.9);
eigs.sym(y, 3, "LM", sigma = 2);
eigs.sym(y, 3, "LM", sigma = 0);
eigs.sym(y, 3, "SM", sigma = 0);
ysy = as(y, "dsyMatrix");
eigs(ysy, 3, "LM", sigma = 1.9);


# TODO: test singularity when using shift-invert
# TODO: eigs("matrix") fails when which == "SM"
# TODO: eigs("matrix") return strange result when which == "SI"

library(rarpack);
x = matrix(rnorm(50), 10);
svd(x);
svds(x, 3);
svd(t(x));
svds(t(x), 3);

