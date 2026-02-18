# Compute Phenotypic Values with Environmental Noise

Adds environmental variation to genetic values to create phenotypic
observations.

\*\*IMPORTANT\*\*: Environmental variance should be calculated ONCE at
generation 0 and held constant across all cycles. Recalculating var_e
each generation artificially maintains heritability and overestimates
genetic gain (fails to model the Bulmer effect correctly).

## Usage

``` r
.compute_phenotypes(genetic_values, environmental_variance)
```

## Arguments

- genetic_values:

  Matrix of genetic values (n_individuals x n_traits)

- environmental_variance:

  Vector of environmental variances (constant across cycles)

## Value

Matrix of phenotypic values (n_individuals x n_traits)
