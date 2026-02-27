# Linear Marker Selection Index (LMSI)

Implements the LMSI which combines phenotypic information with molecular
marker scores from statistically significant markers (Lande & Thompson,
1990). The index is I = b_y' \* y + b_s' \* s, where y are phenotypes
and s are marker scores.

## Usage

``` r
lmsi(
  phen_mat = NULL,
  marker_scores = NULL,
  pmat,
  gmat,
  G_s = NULL,
  wmat,
  wcol = 1,
  selection_intensity = 2.063,
  GAY = NULL
)
```

## Arguments

- phen_mat:

  Matrix of phenotypes (n_genotypes x n_traits). Can be NULL if G_s is
  provided directly (theoretical case where covariance structure is
  known without needing empirical data).

- marker_scores:

  Matrix of marker scores (n_genotypes x n_traits). These are computed
  as s_j = sum(x_jk \* beta_jk) where x_jk is the coded marker value and
  beta_jk is the estimated marker effect for trait j. Can be NULL if G_s
  is provided directly.

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- G_s:

  Genomic covariance matrix explained by markers (n_traits x n_traits).
  This represents Var(s) which approximates Cov(y, s) when markers fully
  explain genetic variance. If provided, phen_mat and marker_scores
  become optional as the covariance structure is specified directly. If
  NULL, computed empirically from marker_scores and phen_mat.

- wmat:

  Economic weights matrix (n_traits x k), or vector.

- wcol:

  Weight column to use if wmat has multiple columns (default: 1).

- selection_intensity:

  Selection intensity k (default: 2.063 for 10% selection).

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation.

## Value

List of class `"lmsi"` with components:

- `b_y`:

  Coefficients for phenotypes (n_traits vector).

- `b_s`:

  Coefficients for marker scores (n_traits vector).

- `b_combined`:

  Combined coefficient vector \[b_y; b_s\] (2\*n_traits vector).

- `P_L`:

  Combined phenotypic-marker covariance matrix (2\*n_traits x
  2\*n_traits).

- `G_L`:

  Combined genetic-marker covariance matrix (2\*n_traits x n_traits).

- `G_s`:

  Genomic covariance matrix explained by markers (n_traits x n_traits).

- `rHI`:

  Index accuracy (correlation between index and breeding objective).

- `sigma_I`:

  Standard deviation of the index.

- `R`:

  Selection response (k \* sigma_I).

- `Delta_H`:

  Expected genetic gain per trait (vector of length n_traits).

- `GA`:

  Overall genetic advance in breeding objective.

- `PRE`:

  Percent relative efficiency (if GAY provided).

- `hI2`:

  Index heritability.

- `summary`:

  Data frame with coefficients and metrics (combined view).

- `phenotype_coeffs`:

  Data frame with phenotype coefficients only.

- `marker_coeffs`:

  Data frame with marker score coefficients only.

- `coeff_analysis`:

  Data frame with coefficient distribution analysis.

## Details

**Mathematical Formulation:**

The LMSI maximizes the correlation between the index \\I\_{LMSI} =
\mathbf{b}\_y^{\prime}\mathbf{y} + \mathbf{b}\_s^{\prime}\mathbf{s}\\
and the breeding objective \\H = \mathbf{w}^{\prime}\mathbf{g}\\.

Combined covariance matrices: \$\$\mathbf{P}\_L = \begin{bmatrix}
\mathbf{P} & \text{Cov}(\mathbf{y}, \mathbf{s}) \\
\text{Cov}(\mathbf{y}, \mathbf{s})^{\prime} & \text{Var}(\mathbf{s})
\end{bmatrix}\$\$ \$\$\mathbf{G}\_L = \begin{bmatrix} \mathbf{G} \\
\mathbf{G}\_s \end{bmatrix}\$\$

where \\\mathbf{P}\\ is the phenotypic variance,
\\\text{Cov}(\mathbf{y}, \mathbf{s})\\ is the covariance between
phenotypes and marker scores (computed from data),
\\\text{Var}(\mathbf{s})\\ is the variance of marker scores,
\\\mathbf{G}\\ is the genotypic variance, and \\\mathbf{G}\_s\\
represents the genetic covariance explained by markers.

Index coefficients: \$\$\mathbf{b}\_{LMSI} = \mathbf{P}\_L^{-1}
\mathbf{G}\_L \mathbf{w}\$\$

Accuracy: \$\$\rho\_{HI} = \sqrt{\frac{\mathbf{b}\_{LMSI}^{\prime}
\mathbf{G}\_L \mathbf{w}}{\mathbf{w}^{\prime} \mathbf{G}
\mathbf{w}}}\$\$

Selection response: \$\$R\_{LMSI} = k \sigma\_{I\_{LMSI}} = k
\sqrt{\mathbf{b}\_{LMSI}^{\prime} \mathbf{P}\_L \mathbf{b}\_{LMSI}}\$\$

Expected genetic gain per trait: \$\$\mathbf{E}\_{LMSI} = k
\frac{\mathbf{G}\_L^{\prime}
\mathbf{b}\_{LMSI}}{\sigma\_{I\_{LMSI}}}\$\$

## References

Lande, R., & Thompson, R. (1990). Efficiency of marker-assisted
selection in the improvement of quantitative traits. Genetics, 124(3),
743-756.

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 4.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load data
data(seldata)
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Simulate marker scores (in practice, computed from QTL mapping)
set.seed(123)
n_genotypes <- 100
n_traits <- ncol(gmat)
marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
  nrow = n_genotypes, ncol = n_traits
)
colnames(marker_scores) <- colnames(gmat)

# Simulate phenotypes
phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
  nrow = n_genotypes, ncol = n_traits
)
colnames(phen_mat) <- colnames(gmat)

# Economic weights
weights <- c(10, 5, 3, 3, 5, 8, 4)

# Calculate LMSI
result <- lmsi(phen_mat, marker_scores, pmat, gmat,
  G_s = NULL, wmat = weights
)
print(result$summary)
} # }
```
