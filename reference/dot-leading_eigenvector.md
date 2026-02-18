# Select the leading eigenvector from a real-eigendecomposition

Returns the eigenvector paired with the \*largest real\* eigenvalue that
is positive (negative eigenvalues indicate numerical noise or rank
deficiency). Normalises the sign so that the first non-zero element is
positive.

## Usage

``` r
.leading_eigenvector(mat, tol = 1e-08)
```

## Arguments

- mat:

  Square matrix to decompose

- tol:

  Eigenvalue tolerance (default 1e-8)

## Value

List: `vector`, `value`, `all_values`
