# Multistage Predetermined Proportional Gain Linear Genomic Selection Index (MPPG-LGSI)

Implements the two-stage Predetermined Proportional Gain LGSI where
breeders specify desired proportional gains between traits at each stage
using GEBVs.

## Usage

``` r
mppg_lgsi(
  Gamma1,
  Gamma,
  A1,
  A,
  C,
  G1,
  P1,
  wmat,
  wcol = 1,
  d1,
  d2,
  U1 = NULL,
  U2 = NULL,
  selection_proportion = 0.1,
  use_young_method = FALSE,
  k1_manual = 2.063,
  k2_manual = 2.063,
  tau = NULL
)
```

## Arguments

- Gamma1:

  GEBV variance-covariance matrix for stage 1 traits (n1 x n1)

- Gamma:

  GEBV variance-covariance matrix for all traits at stage 2 (n x n)

- A1:

  Covariance matrix between GEBVs and true breeding values for stage 1
  (n1 x n1)

- A:

  Covariance matrix between GEBVs and true breeding values for stage 2
  (n x n1)

- C:

  Genotypic variance-covariance matrix for all traits (n x n)

- G1:

  Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)

- P1:

  Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)

- wmat:

  Economic weights vector or matrix (n x k)

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- d1:

  Vector of desired proportional gains for stage 1 (length n1)

- d2:

  Vector of desired proportional gains for stage 2 (length n)

- U1:

  Constraint matrix for stage 1 (n1 x r1), optional

- U2:

  Constraint matrix for stage 2 (n x r2), optional

- selection_proportion:

  Proportion selected at each stage (default: 0.1)

- use_young_method:

  Logical. Use Young's method for selection intensities (default:
  FALSE). Young's method tends to overestimate intensities; manual
  intensities are recommended.

- k1_manual:

  Manual selection intensity for stage 1

- k2_manual:

  Manual selection intensity for stage 2

- tau:

  Standardized truncation point

## Value

List with components similar to mlgsi, plus:

- `beta_P1` - PPG genomic stage 1 coefficients

- `beta_P2` - PPG genomic stage 2 coefficients

- `b_P1` - PPG phenotypic stage 1 coefficients (used for C\* adjustment)

- `theta1` - Proportionality constant for stage 1

- `theta2` - Proportionality constant for stage 2

- `gain_ratios_1` - Achieved gain ratios at stage 1

- `gain_ratios_2` - Achieved gain ratios at stage 2

## Details

**Mathematical Formulation:**

The PPG genomic coefficients are: \$\$\mathbf{\beta}\_{P_1} =
\mathbf{\beta}\_{R_1} + \theta_1
\mathbf{U}\_1(\mathbf{U}\_1'\mathbf{\Gamma}\_1\mathbf{U}\_1)^{-1}\mathbf{d}\_1\$\$
\$\$\mathbf{\beta}\_{P_2} = \mathbf{\beta}\_{R_2} + \theta_2
\mathbf{U}\_2(\mathbf{U}\_2'\mathbf{\Gamma}\mathbf{U}\_2)^{-1}\mathbf{d}\_2\$\$

where proportionality constants are: \$\$\theta_1 =
\frac{\mathbf{d}\_1'(\mathbf{U}\_1'\mathbf{\Gamma}\_1\mathbf{U}\_1)^{-1}\mathbf{U}\_1'\mathbf{A}\_1\mathbf{w}}{\mathbf{d}\_1'(\mathbf{U}\_1'\mathbf{\Gamma}\_1\mathbf{U}\_1)^{-1}\mathbf{d}\_1}\$\$

**Covariance Adjustment:**

The genetic covariance matrix \\\mathbf{C}^\*\\ is adjusted using
phenotypic PPG coefficients \\\mathbf{b}\_{P1} =
\mathbf{P}\_1^{-1}\mathbf{G}\_1\mathbf{P}\_1^{-1}\mathbf{d}\_1\\, which
reflect the same proportional gain constraints as the genomic
coefficients \\\mathbf{\beta}\_{P1}\\. This ensures the adjustment
reflects the actual PPG selection occurring at stage 1.

**Important:** When using custom `U1` matrices (subset constraints), the
phenotypic proxy \\\mathbf{b}\_{P1}\\ uses the standard Tallis formula
(all traits constrained), while the genomic index
\\\mathbf{\beta}\_{P1}\\ respects the `U1` subset. This may cause
\\\mathbf{C}^\*\\ to be slightly over-adjusted. For exact adjustment,
use `U1 = NULL` (default, all traits constrained). Calculating the exact
restricted phenotypic proxy would require implementing the full
MPPG-LPSI projection matrix method for the phenotypic coefficients.

**Note:** Input covariance matrices (`C`, `P1`, `Gamma1`, `Gamma`)
should be positive definite. Non-positive definite matrices may lead to
invalid results or warnings.

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 9,
Section 9.6.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage proportional gain genomic selection
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

reliability <- 0.7
Gamma1 <- reliability * gmat[1:3, 1:3]
Gamma <- reliability * gmat
A1 <- reliability * gmat[1:3, 1:3]
A <- gmat[, 1:3]

# Desired proportional gains
d1 <- c(2, 1, 1)
d2 <- c(3, 2, 1, 1, 1, 0.5, 0.5)

weights <- c(10, 8, 6, 4, 3, 2, 1)

result <- mppg_lgsi(
  Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
  C = gmat, G1 = gmat[1:3, 1:3], P1 = pmat[1:3, 1:3],
  wmat = weights, d1 = d1, d2 = d2
)
} # }
```
