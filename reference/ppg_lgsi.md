# Predetermined Proportional Gains Linear Genomic Selection Index (PPG-LGSI)

Implements the PPG-LGSI where breeders specify desired proportional
gains between traits rather than restricting specific traits to zero.
This is genomic version of PPG-LPSI using GEBVs only.

## Usage

``` r
ppg_lgsi(
  Gamma,
  d,
  wmat = NULL,
  wcol = 1,
  U = NULL,
  k_I = 2.063,
  L_G = 1,
  gmat = NULL,
  GAY = NULL
)
```

## Arguments

- Gamma:

  GEBV variance-covariance matrix (n_traits x n_traits)

- d:

  Vector of desired proportional gains (length n_traits or
  n_constraints). If length n_traits, constraints are applied to all
  traits. If length n_constraints, must provide U matrix.

- wmat:

  Optional. Economic weights for GA/PRE calculation

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- U:

  Optional. Constraint matrix (n_traits x n_constraints). If NULL,
  assumes d applies to all traits (U = I).

- k_I:

  Selection intensity (default: 2.063)

- L_G:

  Standardization constant (default: 1)

- gmat:

  Optional. True genetic variance-covariance matrix for exact accuracy
  calculation. If NULL, uses Gamma as approximation.

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with:

- `summary` - Data frame with coefficients and metrics

- `b` - Vector of PPG-LGSI coefficients (β_PG)

- `E` - Named vector of expected genetic gains per trait

- `theta_G` - Proportionality constant

- `gain_ratios` - Ratios of achieved to desired gains

## Details

**Mathematical Formulation (Chapter 6, Section 6.2):**

Alternative form: \\beta_PG = beta_RG + theta_G \* U \* (U' \* Gamma \*
U)^{-1} \* d\\

Where: - beta_RG = Restricted index coefficients (from RLGSI) - theta_G
= Proportionality constant - d = Vector of desired proportional gains

Proportionality constant: \$\$theta_G = (d' \* (U' \* Gamma \* U)^{-1}
\* U' \* Gamma \* w) / (d' \* (U' \* Gamma \* U)^{-1} \* d)\$\$

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate GEBV variance-covariance matrix
set.seed(123)
n_traits <- 5
Gamma <- matrix(rnorm(n_traits^2), n_traits, n_traits)
Gamma <- (Gamma + t(Gamma)) / 2
diag(Gamma) <- abs(diag(Gamma)) + 2

# Desired proportional gains (e.g., 2:1:1:0:0 ratio)
d <- c(2, 1, 1, 0, 0)

# Economic weights
w <- c(10, 8, 6, 4, 2)

result <- ppg_lgsi(Gamma, d, wmat = w)
print(result$summary)
print(result$gain_ratios) # Should be approximately proportional to d
} # }
```
