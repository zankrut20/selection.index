# Multistage Linear Genomic Selection Index (MLGSI)

Implements the two-stage Linear Genomic Selection Index where selection
is based on GEBVs at both stages with covariance adjustments due to
selection effects.

## Usage

``` r
mlgsi(
  Gamma1,
  Gamma,
  A1,
  A,
  C,
  G1,
  P1,
  wmat,
  wcol = 1,
  selection_proportion = 0.1,
  use_young_method = FALSE,
  k1_manual = 2.063,
  k2_manual = 2.063,
  tau = NULL,
  reliability = NULL
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

- reliability:

  Optional reliability vector for computing A matrices

## Value

List with components:

- `beta1` - Stage 1 genomic index coefficients

- `w` - Economic weights (stage 2 coefficients)

- `stage1_metrics` - List with stage 1 metrics (R1, E1, rho_HI1)

- `stage2_metrics` - List with stage 2 metrics (R2, E2, rho_HI2)

- `Gamma_star` - Adjusted genomic covariance matrix at stage 2

- `C_star` - Adjusted genotypic covariance matrix at stage 2

- `rho_I1I2` - Correlation between stage 1 and stage 2 indices

- `k1` - Selection intensity at stage 1

- `k2` - Selection intensity at stage 2

- `summary_stage1` - Data frame with stage 1 summary

- `summary_stage2` - Data frame with stage 2 summary

## Details

**Mathematical Formulation:**

Stage 1: The genomic index is \\I_1 = \mathbf{\beta}\_1'
\mathbf{\gamma}\_1\\

Coefficients: \\\mathbf{\beta}\_1 =
\mathbf{\Gamma}\_1^{-1}\mathbf{A}\_1\mathbf{w}\_1\\

Stage 2: The index uses economic weights directly: \\I_2 = \mathbf{w}'
\mathbf{\gamma}\\

Adjusted genomic covariance matrix: \$\$\mathbf{\Gamma}^\* =
\mathbf{\Gamma} - u
\frac{\mathbf{A}\_1'\mathbf{\beta}\_1\mathbf{\beta}\_1'\mathbf{A}\_1}{\mathbf{\beta}\_1'\mathbf{\Gamma}\_1\mathbf{\beta}\_1}\$\$

Adjusted genotypic covariance matrix: \$\$\mathbf{C}^\* = \mathbf{C} - u
\frac{\mathbf{G}\_1'\mathbf{b}\_1\mathbf{b}\_1'\mathbf{G}\_1}{\mathbf{b}\_1'\mathbf{P}\_1\mathbf{b}\_1}\$\$

where \\u = k_1(k_1 - \tau)\\

Accuracy at stage 1: \\\rho\_{HI_1} =
\sqrt{\frac{\mathbf{\beta}\_1'\mathbf{\Gamma}\_1\mathbf{\beta}\_1}{\mathbf{w}'\mathbf{C}\mathbf{w}}}\\

Accuracy at stage 2: \\\rho\_{HI_2} =
\sqrt{\frac{\mathbf{w}'\mathbf{\Gamma}^\*\mathbf{w}}{\mathbf{w}'\mathbf{C}^\*\mathbf{w}}}\\

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 9,
Section 9.4.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage genomic selection example
# Stage 1: Select based on GEBVs for 3 traits
# Stage 2: Select based on GEBVs for all 7 traits

# Compute covariance matrices
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Simulate GEBV covariances (in practice, compute from genomic prediction)
set.seed(123)
reliability <- 0.7
Gamma1 <- reliability * gmat[1:3, 1:3]
Gamma <- reliability * gmat
A1 <- reliability * gmat[1:3, 1:3]
A <- gmat[, 1:3]

# Economic weights
weights <- c(10, 8, 6, 4, 3, 2, 1)

# Run MLGSI
result <- mlgsi(
  Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
  C = gmat, G1 = gmat[1:3, 1:3], P1 = pmat[1:3, 1:3],
  wmat = weights, selection_proportion = 0.1
)

print(result$summary_stage1)
print(result$summary_stage2)
} # }
```
