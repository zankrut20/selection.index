# Compute Sum of Products Between Grouped Sums

Efficiently computes sum of products for grouped sum vectors. Equivalent
to crossprod(sums1, sums2) in R.

## Usage

``` r
cpp_crossprod_divided(sums1, sums2, divisor)
```

## Arguments

- sums1:

  Matrix of grouped sums (n_groups x n_traits)

- sums2:

  Matrix of grouped sums (n_groups x n_traits)

- divisor:

  Scalar to divide sums by (e.g., n_replications)

## Value

Matrix of sum of products (n_traits x n_traits)
