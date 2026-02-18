# Multistage Predetermined Proportional Gain Linear Phenotypic Selection Index (MPPG-LPSI)

Implements the two-stage Predetermined Proportional Gain LPSI where
breeders specify desired proportional gains between traits at each
stage.

## Usage

``` r
mppg_lpsi(
  P1,
  P,
  G1,
  C,
  wmat,
  wcol = 1,
  d1,
  d2,
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

- d1:

  Vector of desired proportional gains for stage 1 (length n1)

- d2:

  Vector of desired proportional gains for stage 2 (length n)

- stage1_indices:

  Integer vector specifying which traits correspond to stage 1 (default:
  1:nrow(P1))

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

List with components similar to mlpsi, plus:

- `b_M1` - PPG stage 1 coefficients

- `b_M2` - PPG stage 2 coefficients

- `b_R1` - Restricted stage 1 coefficients

- `b_R2` - Restricted stage 2 coefficients

- `K_M1` - PPG projection matrix for stage 1

- `K_M2` - PPG projection matrix for stage 2

- `theta1` - Proportionality constant for stage 1

- `theta2` - Proportionality constant for stage 2

- `gain_ratios_1` - Achieved gain ratios at stage 1

- `gain_ratios_2` - Achieved gain ratios at stage 2

## Details

**Mathematical Formulation (Chapter 9.3.1, Eq 9.17):**

The PPG coefficients are computed using the projection matrix method:
\$\$\mathbf{b}\_{M_1} = \mathbf{b}\_{R_1} + \theta_1
\mathbf{U}\_1(\mathbf{U}\_1'\mathbf{G}\_1\mathbf{P}\_1^{-1}\mathbf{G}\_1\mathbf{U}\_1)^{-1}\mathbf{d}\_1\$\$
\$\$\mathbf{b}\_{M_2} = \mathbf{b}\_{R_2} + \theta_2
\mathbf{U}\_2(\mathbf{U}\_2'\mathbf{C}\mathbf{P}^{-1}\mathbf{C}\mathbf{U}\_2)^{-1}\mathbf{d}\_2\$\$

where:

- \\\mathbf{b}\_{R_i} = \mathbf{K}\_{M_i}\mathbf{b}\_i\\ are restricted
  coefficients

- \\\mathbf{K}\_{M_i} = \mathbf{I} - \mathbf{Q}\_{M_i}\\ is the
  projection matrix

- \\\theta_i\\ is the proportionality constant computed from
  \\\mathbf{d}\_i\\

- \\\mathbf{U}\_i = \mathbf{I}\\ (all traits constrained)

\$\$\mathbf{b}\_{M_1} = \mathbf{K}\_{M_1} \mathbf{b}\_1\$\$
\$\$\mathbf{b}\_{M_2} = \mathbf{K}\_{M_2} \mathbf{b}\_2\$\$

where \\\mathbf{K}\_{M_i}\\ is computed to achieve proportional gains
specified by \\\mathbf{d}\_i\\

## References

Tallis, G. M. (1962). A selection index for optimum genotype.
Biometrics, 18(1), 120-122.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage proportional gain selection
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

P1 <- pmat[1:3, 1:3]
G1 <- gmat[1:3, 1:3]
P <- pmat
C <- gmat

# Desired proportional gains
d1 <- c(2, 1, 1)  # Trait 1 gains twice as much at stage 1
d2 <- c(3, 2, 1, 1, 1, 0.5, 0.5)  # Different proportions at stage 2

weights <- c(10, 8, 6, 4, 3, 2, 1)

result <- mppg_lpsi(P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
                    d1 = d1, d2 = d2)
} # }
```
