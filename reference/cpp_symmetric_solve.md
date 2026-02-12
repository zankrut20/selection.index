# Solve Symmetric Linear System

Solves Ax = b for symmetric positive definite matrix A using LDLT
decomposition. More efficient than general solve() for symmetric
matrices.

## Usage

``` r
cpp_symmetric_solve(A, b)
```

## Arguments

- A:

  Symmetric positive definite matrix

- b:

  Right-hand side vector

## Value

Solution vector x
