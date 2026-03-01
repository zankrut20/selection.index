# Genetic-Genomic Variance-Covariance Matrix (A)

Computes the genetic-genomic covariance matrix (A) as defined in Chapter
8 (Equation 8.12) for GESIM and related genomic eigen selection indices.

Structure: A = \[\[C, C_g-gamma\], \[C_gamma-g, \\\Gamma\\\]\] (2t x 2t,
square symmetric)

where: - C = Var(g) = true genotypic variance-covariance (t x t) -
\\\Gamma\\ = Var(\\\gamma\\) = genomic variance-covariance (t x t) -
C_g-gamma = Cov(g, \\\gamma\\) = covariance between true BVs and GEBVs
(t x t) - C_gamma-g = Cov(\\\gamma\\, g) = transpose of C_g-gamma (t x
t)

## Usage

``` r
genetic_genomic_varcov(
  gmat,
  Gamma = NULL,
  reliability = NULL,
  C_gebv_g = NULL,
  square = TRUE
)
```

## Arguments

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits)

- Gamma:

  Genomic variance-covariance matrix (n_traits x n_traits). If NULL,
  assumed equal to gmat (perfect prediction).

- reliability:

  Optional. Reliability of GEBVs (r² = squared correlation between GEBV
  and true BV). Can be: - Single value (applied to all traits) - Vector
  of length n_traits (one per trait) - NULL (default): assumes
  C_g,\\\gamma\\ = Gamma (unbiased GEBVs with reliability = 1)

- C_gebv_g:

  Optional. Direct specification of Cov(\\\gamma\\, g) matrix (t x t).
  If provided, overrides reliability parameter.

- square:

  Logical. If TRUE (default), returns (2t × 2t) square matrix as
  required for GESIM. If FALSE, returns (2t × t) rectangular form for
  LMSI.

## Value

Genetic-genomic covariance matrix: - If square = TRUE: (2t × 2t)
symmetric matrix for GESIM/eigen indices - If square = FALSE: (2t × t)
rectangular matrix for LMSI where t is the number of traits

## Details

The genetic-genomic matrix relates selection on phenotypes + GEBVs to
expected genetic gains.

\*\*For GESIM (Chapter 8):\*\* Requires the full (2t × 2t) square matrix
for the eigenproblem: (\\\Phi\\^(-1) A - \\\lambda\\I)b = 0

\*\*For LMSI/CLGSI (Chapter 4):\*\* Can use the rectangular (2t × t)
form in the equation: b = P^(-1) G w, where G is (2t × t).

When reliability is provided: - \\C\_{\gamma g}\\ = diag(\\\sqrt{r^2}\\)

When reliability is NULL: - \\C\_{\gamma g}\\ = Gamma (assumes unbiased
GEBVs, perfect prediction)

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapters 4 &
8.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Simulate genomic covariance
Gamma <- gmat * 0.8

# For GESIM: Get square (2t × 2t) matrix
A_square <- genetic_genomic_varcov(gmat, Gamma, reliability = 0.7)
print(dim(A_square)) # Should be 14 x 14 (2t × 2t)

# For LMSI: Get rectangular (2t × t) matrix
A_rect <- genetic_genomic_varcov(gmat, Gamma, reliability = 0.7, square = FALSE)
print(dim(A_rect)) # Should be 14 x 7 (2t × t)
} # }
```
