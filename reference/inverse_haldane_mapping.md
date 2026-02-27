# Inverse Haldane Mapping Function

Converts recombination fraction back to genetic distance (in Morgans).
This is the inverse of Haldane's mapping function.

## Usage

``` r
inverse_haldane_mapping(recombination_fraction)
```

## Arguments

- recombination_fraction:

  Recombination fraction (r) between 0 and 0.5.

## Value

Genetic distance in Morgans.

## Details

**Mathematical Formula:**

Solving Haldane's equation for d: \$\$d = -\frac{1}{2} \ln(1 - 2r)\$\$

## Examples

``` r
# Convert recombination fraction to distance
inverse_haldane_mapping(0.25) # Returns ~0.347 Morgans
#> [1] 0.3465736
inverse_haldane_mapping(0.5) # Returns Inf (unlinked)
#> Warning: recombination_fraction = 0.5 corresponds to infinite genetic distance
#> [1] Inf
```
