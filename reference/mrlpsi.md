# Multistage Restricted Linear Phenotypic Selection Index (MRLPSI)

Implements the two-stage Restricted Linear Phenotypic Selection Index
where certain traits are constrained to have zero genetic gain at each
stage.

## Usage

``` r
mrlpsi(
  P1,
  P,
  G1,
  C,
  wmat,
  wcol = 1,
  C1,
  C2,
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

- C1:

  Constraint matrix for stage 1 (n1 x r1)

- C2:

  Constraint matrix for stage 2 (n x r2)

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

- `b_R1` - Restricted stage 1 coefficients

- `b_R2` - Restricted stage 2 coefficients

- `K1` - Restriction matrix for stage 1

- `K2` - Restriction matrix for stage 2

## Details

**Mathematical Formulation:**

The restricted coefficients are computed as: \$\$\mathbf{b}\_{R_1} =
\mathbf{K}\_1 \mathbf{b}\_1\$\$ \$\$\mathbf{b}\_{R_2} = \mathbf{K}\_2
\mathbf{b}\_2\$\$

where \\\mathbf{K}\_1 = \mathbf{I}\_1 - \mathbf{Q}\_1\\ and
\\\mathbf{K}\_2 = \mathbf{I}\_2 - \mathbf{Q}\_2\\

and \\\mathbf{Q}\_i =
\mathbf{P}\_i^{-1}\mathbf{G}\_i\mathbf{C}\_i(\mathbf{C}\_i'\mathbf{G}\_i\mathbf{P}\_i^{-1}\mathbf{G}\_i\mathbf{C}\_i)^{-1}\mathbf{C}\_i'\mathbf{G}\_i\\

## References

Kempthorne, O., & Nordskog, A. W. (1959). Restricted selection indices.
Biometrics, 15(1), 10-19.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage restricted selection
# Restrict trait 1 at stage 1, traits 1 and 3 at stage 2

pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

P1 <- pmat[1:3, 1:3]
G1 <- gmat[1:3, 1:3]
P <- pmat
C <- gmat

# Constraint matrices
C1 <- matrix(0, nrow = 3, ncol = 1)
C1[1, 1] <- 1 # Restrict trait 1 at stage 1

C2 <- matrix(0, nrow = 7, ncol = 2)
C2[1, 1] <- 1 # Restrict trait 1 at stage 2
C2[3, 2] <- 1 # Restrict trait 3 at stage 2

weights <- c(10, 8, 6, 4, 3, 2, 1)

result <- mrlpsi(
  P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
  C1 = C1, C2 = C2
)
} # }
```
