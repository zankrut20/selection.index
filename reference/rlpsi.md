# Restricted Linear Phenotypic Selection Index (RLPSI)

Restricted Linear Phenotypic Selection Index (RLPSI)

## Usage

``` r
rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = NULL, C = NULL, GAY)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix

- gmat:

  Genotypic variance-covariance matrix

- wmat:

  Weight matrix

- wcol:

  Weight column number (default: 1)

- restricted_traits:

  Vector of trait indices to restrict (default: NULL). If provided, a
  constraint matrix C is auto-generated to enforce zero gain on these
  traits. Example: c(1, 3) restricts traits 1 and 3 to zero gain.

- C:

  Constraint matrix (n_traits x n_constraints). Each column is a
  restriction. Alternative to restricted_traits for custom constraints.
  Ignored if restricted_traits is provided.

- GAY:

  Genetic advance of comparative trait (optional)

## Value

List with:

- `summary` - Data frame with coefficients (b.\*), GA, PRE, Delta_G,
  rHI, hI2

- `b` - Numeric vector of selection index coefficients

- `Delta_G` - Named vector of realized correlated responses per trait

- `C` - Constraint matrix used

## Examples

``` r
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
wmat <- weight_mat(weight)

# Easy way: Restrict traits 1 and 3 to zero gain
result <- rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = c(1, 3))

# Advanced way: Provide custom constraint matrix
C <- diag(ncol(pmat))[, 1, drop = FALSE]
result <- rlpsi(pmat, gmat, wmat, wcol = 1, C = C)
```
