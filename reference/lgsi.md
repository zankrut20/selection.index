# Linear Genomic Selection Index (LGSI)

Implements the Linear Genomic Selection Index where selection is based
solely on Genomic Estimated Breeding Values (GEBVs). This is used for
selecting candidates that have been genotyped but not phenotyped (e.g.,
in a testing population).

## Usage

``` r
lgsi(
  gebv_mat,
  gmat,
  wmat,
  wcol = 1,
  reliability = NULL,
  selection_intensity = 2.063,
  GAY = NULL
)
```

## Arguments

- gebv_mat:

  Matrix of GEBVs (n_genotypes x n_traits)

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits)

- wmat:

  Economic weights matrix (n_traits x k), or vector

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- reliability:

  Optional. Reliability of GEBVs (correlation between GEBV and true BV).
  Can be: - Single value (applied to all traits) - Vector of length
  n_traits (one per trait) - NULL (default): estimated from GEBV
  variance (assumes reliability = GEBV_var / G_var)

- selection_intensity:

  Selection intensity i (default: 2.063 for 10% selection)

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with components:

- `b` - Index coefficients

- `P_gebv` - GEBV variance-covariance matrix

- `reliability` - Reliability values used

- `Delta_H` - Expected genetic advance per trait

- `GA` - Overall genetic advance in the index

- `PRE` - Percent relative efficiency (if GAY provided)

- `hI2` - Index heritability

- `rHI` - Index accuracy

- `sigma_I` - Standard deviation of the index

- `summary` - Data frame with all metrics

## Details

**Mathematical Formulation:**

The LGSI maximizes the correlation between the index I = b' \* gebv and

Index coefficients: \\\mathbf{b} = \mathbf{P}\_{\hat{g}}^{-1}
\mathbf{C}\_{\hat{g}g} \mathbf{w}\\

Where: - \\\mathbf{P}\_{\hat{g}}\\ = Var(gebv) - variance-covariance of
GEBVs - \\\mathbf{C}\_{\hat{g}g}\\ = Cov(gebv, g) - covariance between
GEBVs and true breeding values

If reliability (r) is known: \\\mathbf{C}\_{\hat{g}g} = \text{diag}(r)
\mathbf{P}\_{\hat{g}}\\

Expected response: \\\Delta \mathbf{H} = \frac{i}{\sigma_I}
\mathbf{C}\_{\hat{g}g} \mathbf{b}\\

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Simulate GEBVs (in practice, these come from genomic prediction)
set.seed(123)
n_genotypes <- 100
n_traits <- ncol(gmat)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
  nrow = n_genotypes, ncol = n_traits
)
colnames(gebv_mat) <- colnames(gmat)

# Economic weights
weights <- c(10, 5, 3, 3, 5, 8, 4)

# Calculate LGSI
result <- lgsi(gebv_mat, gmat, weights, reliability = 0.7)
print(result$summary)
} # }
```
