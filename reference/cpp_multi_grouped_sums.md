# Compute Multiple Grouped Sums at Once

Computes grouped sums for multiple grouping variables simultaneously.
More efficient than calling cpp_grouped_sums multiple times.

## Usage

``` r
cpp_multi_grouped_sums(data_mat, group_indices)
```

## Arguments

- data_mat:

  Numeric matrix (n_obs x n_traits)

- group_indices:

  List of integer vectors, each representing a grouping variable

## Value

List of matrices, one for each grouping variable
