# Extract Symmetric Submatrix

Extracts a symmetric submatrix given row/column indices. Used for
selecting trait combinations from full covariance matrices.

## Usage

``` r
cpp_extract_submatrix(mat, indices)
```

## Arguments

- mat:

  Numeric matrix (symmetric)

- indices:

  Integer vector of indices (1-based, converted to 0-based internally)

## Value

Symmetric submatrix
