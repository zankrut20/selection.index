# Combinatorial Linear Phenotypic Selection Index

Build all possible Smith-Hazel selection indices from trait
combinations, with optional exclusion of specific traits.

This function systematically evaluates indices for all combinations of
ncomb traits, which is useful for identifying the most efficient subset
of traits for selection.

## Usage

``` r
lpsi(ncomb, pmat, gmat, wmat, wcol = 1, GAY, excluding_trait = NULL)
```

## Arguments

- ncomb:

  Number of traits per combination

- pmat:

  Phenotypic variance-covariance matrix

- gmat:

  Genotypic variance-covariance matrix

- wmat:

  Weight matrix

- wcol:

  Weight column number if more than one weight set (default: 1)

- GAY:

  Genetic advance of comparative trait (optional)

- excluding_trait:

  Optional. Traits to exclude from combinations. Can be: (1) numeric
  vector of trait indices (e.g., c(1, 3)), (2) character vector of trait
  names (e.g., c("sypp", "dtf")), (3) data frame/matrix columns with
  trait data (trait names extracted from column names). When specified,
  only combinations that do NOT contain any of these traits are
  returned.

## Value

Data frame of all possible selection indices with metrics (GA, PRE,
Delta_G, rHI, hI2)

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
wmat <- weight_mat(weight)

# Build all 3-trait indices
result <- lpsi(ncomb = 3, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1)

# Exclude specific traits
result <- lpsi(
  ncomb = 3, pmat = pmat, gmat = gmat, wmat = wmat,
  excluding_trait = c(1, 3)
)
} # }
```
