# Multistage Linear Phenotypic Selection Index (MLPSI)

Implements the two-stage Linear Phenotypic Selection Index where
selection occurs at two different stages with covariance adjustments due
to selection effects using Cochran/Cunningham's method.

## Usage

``` r
mlpsi(
  P1,
  P,
  G1,
  C,
  wmat,
  wcol = 1,
  stage1_indices = NULL,
  selection_proportion = 0.1,
  use_young_method = FALSE,
  k1_manual = 2.063,
  k2_manual = 2.063,
  tau = NULL
)
```

## Arguments

- P1:

  Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)

- P:

  Phenotypic variance-covariance matrix for all traits at stage 2 (n x
  n)

- G1:

  Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)

- C:

  Genotypic variance-covariance matrix for all traits (n x n)

- wmat:

  Economic weights vector or matrix (n x k)

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- stage1_indices:

  Integer vector specifying which traits (columns of P and C) correspond
  to stage 1. Default is 1:nrow(P1), assuming stage 1 traits are the
  first n1 traits. Use this to specify non-contiguous traits, e.g., c(1,
  3, 5) for traits 1, 3, and 5.

- selection_proportion:

  Proportion selected at each stage (default: 0.1)

- use_young_method:

  Logical. Use Young's method for selection intensities (default:
  FALSE). Young's method tends to overestimate intensities; manual
  intensities are recommended.

- k1_manual:

  Manual selection intensity for stage 1 (used if use_young_method =
  FALSE)

- k2_manual:

  Manual selection intensity for stage 2 (used if use_young_method =
  FALSE)

- tau:

  Standardized truncation point (default: computed from
  selection_proportion)

## Value

List with components:

- `b1` - Stage 1 index coefficients

- `b2` - Stage 2 index coefficients

- `stage1_metrics` - List with stage 1 metrics (R1, E1, rho_H1)

- `stage2_metrics` - List with stage 2 metrics (R2, E2, rho_H2)

- `P_star` - Adjusted phenotypic covariance matrix at stage 2

- `C_star` - Adjusted genotypic covariance matrix at stage 2

- `rho_12` - Correlation between stage 1 and stage 2 indices

- `k1` - Selection intensity at stage 1

- `k2` - Selection intensity at stage 2

- `summary_stage1` - Data frame with stage 1 summary

- `summary_stage2` - Data frame with stage 2 summary

## Details

**Mathematical Formulation:**

Stage 1 index coefficients: \$\$\mathbf{b}\_1 =
\mathbf{P}\_1^{-1}\mathbf{G}\_1\mathbf{w}\$\$

Stage 2 index coefficients: \$\$\mathbf{b}\_2 =
\mathbf{P}^{-1}\mathbf{G}\mathbf{w}\$\$

Adjusted phenotypic covariance matrix (Cochran/Cunningham):
\$\$\mathbf{P}^\* = \mathbf{P} - u
\frac{Cov(\mathbf{y},\mathbf{x}\_1)\mathbf{b}\_1\mathbf{b}\_1'Cov(\mathbf{x}\_1,\mathbf{y})}{\mathbf{b}\_1'\mathbf{P}\_1\mathbf{b}\_1}\$\$

Adjusted genotypic covariance matrix: \$\$\mathbf{C}^\* = \mathbf{C} - u
\frac{\mathbf{G}\_1'\mathbf{b}\_1\mathbf{b}\_1'\mathbf{G}\_1}{\mathbf{b}\_1'\mathbf{P}\_1\mathbf{b}\_1}\$\$

where \\u = k_1(k_1 - \tau)\\

Selection response: \\R_1 = k_1
\sqrt{\mathbf{b}\_1'\mathbf{P}\_1\mathbf{b}\_1}\\, \\R_2 = k_2
\sqrt{\mathbf{b}\_2'\mathbf{P}^\*\mathbf{b}\_2}\\

## References

Cochran, W. G. (1951). Improvement by means of selection. Proceedings of
the Second Berkeley Symposium on Mathematical Statistics and
Probability, 449-470.

Cunningham, E. P. (1975). Multi-stage index selection. Theoretical and
Applied Genetics, 46(2), 55-61.

Young, S. S. Y. (1964). Multi-stage selection for genetic gain.
Heredity, 19(1), 131-144.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage selection example
# Stage 1: Select based on 3 traits
# Stage 2: Select based on all 7 traits

# Compute variance-covariance matrices
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Stage 1 uses first 3 traits
P1 <- pmat[1:3, 1:3]
G1 <- gmat[1:3, 1:3]

# Stage 2 uses all 7 traits
P <- pmat
C <- gmat

# Economic weights
weights <- c(10, 8, 6, 4, 3, 2, 1)

# Run MLPSI (default: stage1_indices = 1:3)
result <- mlpsi(
  P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
  selection_proportion = 0.1
)

# Or with non-contiguous traits (e.g., traits 1, 3, 5 at stage 1):
# P1 <- pmat[c(1,3,5), c(1,3,5)]
# G1 <- gmat[c(1,3,5), c(1,3,5)]
# result <- mlpsi(P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
#                 stage1_indices = c(1, 3, 5))

print(result$summary_stage1)
print(result$summary_stage2)
} # }
```
