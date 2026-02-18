# Adjust genotypic covariance matrix after stage 1 selection

Adjust genotypic covariance matrix after stage 1 selection

## Usage

``` r
.adjust_genotypic_matrix(C, G1, b1, k1, tau, P1, stage1_indices = NULL)
```

## Arguments

- stage1_indices:

  Indices of stage 1 traits in the full matrix (default: 1:nrow(G1))
