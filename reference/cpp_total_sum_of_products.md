# Total Sum of Products

Computes the total sum of products matrix corrected for the mean.
TSP(i,j) = sum(x_i \* x_j) - CF(i,j)

## Usage

``` r
cpp_total_sum_of_products(data_mat, CF)
```

## Arguments

- data_mat:

  Data matrix (n_obs x n_traits)

- CF:

  Correction factor matrix

## Value

Symmetric sum of products matrix
