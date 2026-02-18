# Adjust phenotypic covariance matrix after stage 1 selection

Adjust phenotypic covariance matrix after stage 1 selection

## Usage

``` r
.adjust_phenotypic_matrix(P, P1, b1, k1, tau, stage1_indices = NULL)
```

## Arguments

- stage1_indices:

  Indices of stage 1 traits in the full matrix (default: 1:nrow(P1))
