# Compute Grouped Sums for Matrix Columns

Efficiently computes grouped sums for all columns of a matrix.
Equivalent to rowsum() in R but optimized for multiple columns.

## Usage

``` r
cpp_grouped_sums(data_mat, group_idx)
```

## Arguments

- data_mat:

  Numeric matrix (n_obs x n_traits)

- group_idx:

  Integer vector of group indices (1-based, converted to 0-based
  internally)

## Value

Matrix of grouped sums (n_groups x n_traits)
