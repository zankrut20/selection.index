# Compute Correction Factor Matrix

Computes correction factor for all trait pairs. CF\[i,j\] = (sum_i \*
sum_j) / n_obs

## Usage

``` r
cpp_correction_factor_matrix(data_mat)
```

## Arguments

- data_mat:

  Numeric matrix (n_obs x n_traits)

## Value

Matrix of correction factors (n_traits x n_traits)
