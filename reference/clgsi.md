# Combined Linear Genomic Selection Index (CLGSI)

Implements the Combined Linear Genomic Selection Index where selection
combines both phenotypic observations and Genomic Estimated Breeding
Values (GEBVs). This is used for selecting candidates with both
phenotype and genotype data (e.g., in a training population).

## Usage

``` r
clgsi(
  phen_mat = NULL,
  gebv_mat = NULL,
  pmat,
  gmat,
  wmat,
  wcol = 1,
  P_y = NULL,
  P_g = NULL,
  P_yg = NULL,
  reliability = NULL,
  selection_intensity = 2.063,
  GAY = NULL
)
```

## Arguments

- phen_mat:

  Matrix of phenotypes (n_genotypes x n_traits). Can be NULL if P_y and
  P_yg are provided.

- gebv_mat:

  Matrix of GEBVs (n_genotypes x n_traits). Can be NULL if P_g and P_yg
  are provided.

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits)

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits)

- wmat:

  Economic weights matrix (n_traits x k), or vector

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- P_y:

  Optional. Phenotypic variance-covariance matrix computed from data
  (n_traits x n_traits). If NULL (default), uses pmat or computes from
  phen_mat using cov().

- P_g:

  Optional. GEBV variance-covariance matrix (n_traits x n_traits). If
  NULL (default), computed from gebv_mat using cov().

- P_yg:

  Optional. Covariance matrix between phenotypes and GEBVs (n_traits x
  n_traits). If NULL (default), computed from phen_mat and gebv_mat
  using cov().

- reliability:

  Optional. Reliability of GEBVs (r^2, the squared correlation). See
  lgsi() for details.

- selection_intensity:

  Selection intensity i (default: 2.063 for 10% selection)

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with components:

- `b_y` - Coefficients for phenotypes

- `b_g` - Coefficients for GEBVs

- `b_combined` - Full coefficient vector \[b_y; b_g\]

- `P_combined` - Combined variance matrix

- `Delta_H` - Expected genetic advance per trait

- `GA` - Overall genetic advance

- `PRE` - Percent relative efficiency

- `hI2` - Index heritability

- `rHI` - Index accuracy

- `summary` - Data frame with all metrics

## Details

**Mathematical Formulation:**

The CLGSI combines phenotypic and genomic information: \$\$I =
\mathbf{b}\_y' \mathbf{y} + \mathbf{b}\_g' \hat{\mathbf{g}}\$\$

Coefficients are obtained by solving the partitioned system:
\$\$\begin{bmatrix} \mathbf{b}\_y \\ \mathbf{b}\_g \end{bmatrix} =
\begin{bmatrix} \mathbf{P} & \mathbf{P}\_{y\hat{g}} \\
\mathbf{P}\_{y\hat{g}}' & \mathbf{P}\_{\hat{g}} \end{bmatrix}^{-1}
\begin{bmatrix} \mathbf{G} \\ \mathbf{C}\_{\hat{g}g} \end{bmatrix}
\mathbf{w}\$\$

Where: - \\\mathbf{P}\\ = Var(phenotypes) - \\\mathbf{P}\_{\hat{g}}\\ =
Var(GEBVs) - \\\mathbf{P}\_{y\hat{g}}\\ = Cov(phenotypes, GEBVs) -
\\\mathbf{G}\\ = Genotypic variance-covariance -
\\\mathbf{C}\_{\hat{g}g}\\ = Cov(GEBV, true BV)

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate phenotypes and GEBVs
set.seed(123)
n_genotypes <- 100
n_traits <- ncol(gmat)

phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                   nrow = n_genotypes, ncol = n_traits)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
                   nrow = n_genotypes, ncol = n_traits)
colnames(phen_mat) <- colnames(gmat)
colnames(gebv_mat) <- colnames(gmat)

# Economic weights
weights <- c(10, 5, 3, 3, 5, 8, 4)

# Calculate CLGSI
result <- clgsi(phen_mat, gebv_mat, pmat, gmat, weights, reliability = 0.7)
print(result$summary)
} # }
```
