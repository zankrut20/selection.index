# Grouped Sum of Products

Computes the sum of products for grouped data. GSP(i,j) = sum_g
\[(sum_i_g \* sum_j_g) / n_g\] - CF(i,j)

## Usage

``` r
cpp_grouped_sum_of_products(group_sums, group_counts, CF)
```

## Arguments

- group_sums:

  Matrix of group sums (n_groups x n_traits)

- group_counts:

  Vector of group sizes

- CF:

  Correction factor matrix

## Value

Symmetric sum of products matrix
