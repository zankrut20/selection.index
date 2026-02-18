# Compute standard metrics for an eigen-based index

Computes all standard metrics using C++ math primitives. For ESIM family
indices the eigenvalue lambda^2 IS the index heritability.

## Usage

``` r
.eigen_index_metrics(b, P, G, lambda2 = NULL, k_I = 2.063)
```

## Arguments

- b:

  Index coefficient vector (eigenvector or transformed eigenvector)

- P:

  Phenotypic VCV matrix

- G:

  Genotypic VCV matrix (C in LaTeX notation)

- lambda2:

  Eigenvalue associated with b (used for h^2_I when exact)

- k_I:

  Selection intensity constant (default 2.063)
