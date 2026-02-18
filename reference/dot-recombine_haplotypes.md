# Recombine Two Haplotypes

Simulates meiotic recombination between two parental haplotypes using
the recombination fractions derived from Haldane's mapping function.

## Usage

``` r
.recombine_haplotypes(haplotype1, haplotype2, recombination_fractions)
```

## Arguments

- haplotype1:

  First parental haplotype (vector of length n_loci)

- haplotype2:

  Second parental haplotype (vector of length n_loci)

- recombination_fractions:

  Vector of recombination fractions between adjacent loci (length
  n_loci - 1)

## Value

Recombinant haplotype (vector of length n_loci)
