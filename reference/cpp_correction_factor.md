# Correction Factor Matrix

Computes the correction factor matrix for ANOVA calculations. CF(i,j) =
(sum_i \* sum_j) / n

## Usage

``` r
cpp_correction_factor(total_sums, n_obs)
```

## Arguments

- total_sums:

  Vector of column sums

- n_obs:

  Number of observations

## Value

Symmetric correction factor matrix
