# Predetermined Proportional Gains (PPG-LPSI)

Implements the PPG-LPSI where breeders specify desired proportional
gains between traits rather than restricting specific traits to zero.
Based on Tallis (1962).

## Usage

``` r
ppg_lpsi(pmat, gmat, k, wmat = NULL, wcol = 1, GAY)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits)

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits)

- k:

  Vector of desired proportional gains (length n_traits). Example: k =
  c(2, 1, 1) means trait 1 should gain twice as much as traits 2 and 3.

- wmat:

  Optional weight matrix for GA/PRE calculation

- wcol:

  Weight column number (default: 1)

- GAY:

  Genetic advance of comparative trait (optional)

## Value

List with:

- `summary` - Data frame with coefficients and metrics

- `b` - Vector of PPG-LPSI coefficients

- `Delta_G` - Expected genetic gains per trait

- `phi` - Proportionality constant

## Details

**Mathematical Formulation (Chapter 3, Section 3.2):**

The PPG-LPSI achieves gains in specific proportions: Delta_G = phi\*k

Coefficient formula (Tallis, 1962): \$\$b = P^{-1}G(G'P^{-1}G)^{-1}k\$\$

Where: - k = Vector of desired proportions - phi = Proportionality
constant (determined by selection intensity and variances)

The constraint ensures Delta_G1:Delta_G2:Delta_G3 = k1:k2:k3

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Gains in ratio 2:1:1:1:1:1:1
k <- c(2, 1, 1, 1, 1, 1, 1)
result <- ppg_lpsi(pmat, gmat, k)
} # }
```
