## R wrapper of ARPACK for large scale eigen value/vector problems

### Introduction

**rARPACK** is an R wrapper of the
[ARPACK library](http://www.caam.rice.edu/software/ARPACK/)
to solve large scale eigen
value/vector problems. It is typically used to compute a few eigen
values/vectors of an `n` by `n` matrix, e.g., the `k` largest eigen values, which
is usually more efficient than `eigen()` if `k << n`. 

Currently this package provides function `eigs()` for eigenvalue/eigenvector
problems, and `svds()` for truncated SVD. Different matrix types in R,
including sparse matrices, are supported. Below is a list of implemented ones:

- `matrix` (defined in base R)
- `dgeMatrix` (defined in **Matrix** package, for general matrices)
- `dsyMatrix` (defined in **Matrix** package, for symmetric matrices)
- `dgCMatrix` (defined in **Matrix** package, for column oriented sparse matrices)

### Example

Set up matrices.

```
library(rARPACK)
library(Matrix)
n = 10
k = 5
set.seed(123)
x = matrix(rnorm(n^2), n)
x[sample(n^2, floor(n^2 / 2))] = 0

xsp = as(x, "sparseMatrix");

y = crossprod(x);
ysy = as(y, "symmetricMatrix");
```

Compute the largest 5 eigenvalues with corresponding eigenvectors.

```
eigs(y, k)
eigs(ysy, k)
```

Compute the smallest 5 eigenvalues with corresponding eigenvectors

```
eigs(ysy, k, which = "SM")
```

or

```
eigs(ysy, k, sigma = 0)
```

Complex eigenvalues, on sparse matrix.

```
eigs(xsp, k)
```

SVD, you can specify the number of singular values(`k`),
number of left singular vectors(`nu`) and number of right singular
vectors(`nv`).

```
n = 10
p = 8
k = 5
set.seed(123)
x = matrix(rnorm(n * p), n)
x[sample(n * p, floor(n * p / 2))] = 0
xsp = as(x, "sparseMatrix")
xspt = t(xsp)

svds(x, k)
svds(xsp, k, nu = 0)
svds(xspt, k, nu = 0, nv = 0)
```

