### R wrapper of ARPACK for large scale eigen value/vector problems

**rarpack** is an R wrapper of the
[ARPACK library](http://www.caam.rice.edu/software/ARPACK/)
to solve large scale eigen
value/vector problems. It is typically used to compute a few eigen
values/vectors of an `n` by `n` matrix, e.g., the `k` largest eigen values, which
is usually more efficient than `eigen()` if `k << n`. Currently this package
provides the function `eigs()` which does the similar job as in Matlab and
Octave.
