# Haldane's Mapping Function

Converts genetic distance (in Morgans) to recombination fraction using
Haldane's mapping function. This function models the relationship
between genetic distance and the probability of recombination between
loci.

## Usage

``` r
haldane_mapping(distance)
```

## Arguments

- distance:

  Genetic distance in Morgans (scalar or vector). One Morgan corresponds
  to a 50% recombination frequency.

## Value

Recombination fraction (r) ranging from 0 to 0.5. - r = 0 indicates
complete linkage (no recombination) - r = 0.5 indicates independent
assortment (unlinked loci)

## Details

**Mathematical Formula (Chapter 10, Section 10.1):**

The relationship between recombination fraction (r) and genetic distance
(d): \$\$r = \frac{1}{2}(1 - e^{-2d})\$\$

Where: - d = Genetic distance in Morgans - r = Recombination fraction
(probability of recombination per meiosis)

This function assumes no crossover interference beyond that implied by
the mapping function itself.

## Examples

``` r
# Zero distance means complete linkage (no recombination)
haldane_mapping(0)  # Returns 0
#> [1] 0

# 1 Morgan distance
haldane_mapping(1)  # Returns ~0.43
#> [1] 0.4323324

# Large distance approaches 0.5 (independent assortment)
haldane_mapping(10)  # Returns ~0.5
#> [1] 0.5

# Vector of distances
distances <- c(0, 0.1, 0.5, 1.0, 2.0)
haldane_mapping(distances)
#> [1] 0.00000000 0.09063462 0.31606028 0.43233236 0.49084218
```
