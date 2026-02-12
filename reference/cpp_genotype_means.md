# Compute Genotype Means Matrix

Efficiently computes means for each genotype across all traits.
Equivalent to rowsum(data, genotypes) / counts but optimized.

## Usage

``` r
cpp_genotype_means(data_mat, gen_idx)
```

## Arguments

- data_mat:

  Numeric matrix (n_obs x n_traits)

- gen_idx:

  Integer vector of genotype indices (1-based)

## Value

Matrix of genotype means (n_genotypes x n_traits)
