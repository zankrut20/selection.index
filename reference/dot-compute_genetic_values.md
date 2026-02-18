# Compute Genetic Values for Population

Calculates breeding values (genetic values) for all individuals in the
population based on their diploid genotypes and QTL effects.

## Usage

``` r
.compute_genetic_values(haplotype1, haplotype2, qtl_effects)
```

## Arguments

- haplotype1:

  First haplotype matrix (n_individuals x n_loci)

- haplotype2:

  Second haplotype matrix (n_individuals x n_loci)

- qtl_effects:

  QTL effect matrix (n_loci x n_traits)

## Value

Matrix of genetic values (n_individuals x n_traits)
