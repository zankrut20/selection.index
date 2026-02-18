# Compute standard metrics for a genomic eigen-based index

Computes all standard metrics using C++ math primitives. For genomic
eigen indices, the eigenvalue lambda^2 IS the index heritability.

## Usage

``` r
.genomic_eigen_index_metrics(b, Phi, A, lambda2 = NULL, k_I = 2.063)
```

## Arguments

- b:

  Index coefficient vector (eigenvector or transformed eigenvector)

- Phi:

  Combined phenotypic/marker variance-covariance matrix

- A:

  Combined genotypic/marker variance-covariance matrix

- lambda2:

  Eigenvalue associated with b (used for h^2_I when exact)

- k_I:

  Selection intensity constant (default 2.063)
