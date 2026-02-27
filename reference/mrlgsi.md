# Multistage Restricted Linear Genomic Selection Index (MRLGSI)

Implements the two-stage Restricted Linear Genomic Selection Index where
certain traits are constrained to have zero genetic gain at each stage
using GEBVs.

## Usage

``` r
mrlgsi(
  Gamma1,
  Gamma,
  A1,
  A,
  C,
  G1,
  P1,
  wmat,
  wcol = 1,
  C1,
  C2,
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

- C1:

  Constraint matrix for stage 1 (n1 x r1)

- C2:

  Constraint matrix for stage 2 (n x r2)

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

- `beta_R1` - Restricted stage 1 coefficients

- `beta_R2` - Restricted stage 2 coefficients

- `K_G1` - Restriction matrix for stage 1

- `K_G2` - Restriction matrix for stage 2

## Details

**Mathematical Formulation:**

The restricted genomic coefficients are: \$\$\mathbf{\beta}\_{R_1} =
\mathbf{K}\_{G_1}\mathbf{\beta}\_1\$\$ \$\$\mathbf{\beta}\_{R_2} =
\mathbf{K}\_{G_2}\mathbf{w}\$\$

where restriction matrices are computed similarly to RLGSI

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 9,
Section 9.5.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage restricted genomic selection
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

reliability <- 0.7
Gamma1 <- reliability * gmat[1:3, 1:3]
Gamma <- reliability * gmat
A1 <- reliability * gmat[1:3, 1:3]
A <- gmat[, 1:3]

# Constraint matrices
C1 <- matrix(0, nrow = 3, ncol = 1)
C1[1, 1] <- 1 # Restrict trait 1 at stage 1

C2 <- matrix(0, nrow = 7, ncol = 2)
C2[1, 1] <- 1 # Restrict trait 1 at stage 2
C2[3, 2] <- 1 # Restrict trait 3 at stage 2

weights <- c(10, 8, 6, 4, 3, 2, 1)

result <- mrlgsi(
  Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
  C = gmat, G1 = gmat[1:3, 1:3], P1 = pmat[1:3, 1:3],
  wmat = weights, C1 = C1, C2 = C2
)
} # }
```
