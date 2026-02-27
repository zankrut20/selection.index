# Genome-Wide Linear Marker Selection Index (GW-LMSI)

Implements the GW-LMSI which uses all available genome-wide markers
directly as predictors in the selection index. Unlike LMSI which uses
aggregated marker scores per trait, GW-LMSI treats each individual
marker as a separate predictor.

## Usage

``` r
gw_lmsi(
  marker_mat,
  trait_mat = NULL,
  gmat,
  P_GW = NULL,
  G_GW = NULL,
  wmat,
  wcol = 1,
  lambda = 0,
  selection_intensity = 2.063,
  GAY = NULL
)
```

## Arguments

- marker_mat:

  Matrix of marker genotypes (n_genotypes x n_markers). Typically coded
  as -1, 0, 1 or 0, 1, 2.

- trait_mat:

  Matrix of trait values (n_genotypes x n_traits). Used to compute G_GW
  if not provided.

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- P_GW:

  Marker covariance matrix (n_markers x n_markers). If NULL, computed as
  Var(marker_mat).

- G_GW:

  Covariance between markers and traits (n_markers x n_traits). If NULL,
  computed as Cov(marker_mat, trait_mat).

- wmat:

  Economic weights matrix (n_traits x k), or vector.

- wcol:

  Weight column to use if wmat has multiple columns (default: 1).

- lambda:

  Ridge regularization parameter (default: 0). If lambda \> 0, uses
  P_GW + lambda\*I for regularization. Automatic warnings issued when
  n_markers \> n_genotypes (high-dimensional case) or when P_GW is
  ill-conditioned. Recommended values: 0.01-0.1 times mean(diag(P_GW)).

- selection_intensity:

  Selection intensity k (default: 2.063 for 10% selection).

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation.

## Value

List of class `"gw_lmsi"` with components:

- `b`:

  Index coefficients for markers (n_markers).

- `P_GW`:

  Marker covariance matrix (n_markers x n_markers).

- `G_GW`:

  Covariance between markers and traits (n_markers x n_traits).

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

- `lambda`:

  Ridge regularization parameter used.

- `n_markers`:

  Number of markers.

- `high_dimensional`:

  Logical indicating if n_markers \> n_genotypes.

- `condition_number`:

  Condition number of P_GW (if computed).

- `summary`:

  Data frame with metrics.

## Details

**Mathematical Formulation:**

The GW-LMSI maximizes the correlation between the index \\I\_{GW} =
\mathbf{b}\_{GW}^{\prime}\mathbf{m}\\ and the breeding objective \\H =
\mathbf{w}^{\prime}\mathbf{g}\\.

Marker covariance matrix: \$\$\mathbf{P}\_{GW} = Var(\mathbf{m})\$\$

Covariance between markers and traits: \$\$\mathbf{G}\_{GW} =
Cov(\mathbf{m}, \mathbf{g})\$\$

Index coefficients (with optional Ridge regularization):
\$\$\mathbf{b}\_{GW} = (\mathbf{P}\_{GW} + \lambda \mathbf{I})^{-1}
\mathbf{G}\_{GW} \mathbf{w}\$\$

Accuracy: \$\$\rho\_{HI} = \sqrt{\frac{\mathbf{b}\_{GW}^{\prime}
\mathbf{G}\_{GW} \mathbf{w}}{\mathbf{w}^{\prime} \mathbf{G}
\mathbf{w}}}\$\$

Selection response: \$\$R\_{GW} = k \sqrt{\mathbf{b}\_{GW}^{\prime}
\mathbf{P}\_{GW} \mathbf{b}\_{GW}}\$\$

Expected genetic gain per trait: \$\$\mathbf{E}\_{GW} = k
\frac{\mathbf{G}\_{GW}^{\prime} \mathbf{b}\_{GW}}{\sigma\_{I\_{GW}}}\$\$

**Note on Singularity Detection and Regularization:** The function
automatically detects problematic cases:

1\. \*\*High-dimensional case\*\*: When n_markers \> n_genotypes, P_GW
is mathematically singular (rank-deficient). The function issues a
warning and suggests an appropriate lambda value.

2\. \*\*Ill-conditioned case\*\*: When P_GW has a high condition number
(\> 1e10), indicating numerical instability.

3\. \*\*Numerical singularity\*\*: When P_GW has eigenvalues near zero.

Ridge regularization adds λI to P_GW, ensuring positive definiteness.
Recommended lambda values are 0.01-0.1 times the average diagonal
element of P_GW. Users can also set lambda = 0 to force generalized
inverse (less stable but sometimes needed).

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
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Simulate marker data
set.seed(123)
n_genotypes <- 100
n_markers <- 200
n_traits <- ncol(gmat)

# Marker matrix (coded as 0, 1, 2)
marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
  nrow = n_genotypes, ncol = n_markers
)

# Trait matrix
trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
  nrow = n_genotypes, ncol = n_traits
)

# Economic weights
weights <- c(10, 5, 3, 3, 5, 8, 4)

# Calculate GW-LMSI with Ridge regularization
result <- gw_lmsi(marker_mat, trait_mat, gmat,
  wmat = weights, lambda = 0.01
)
print(result$summary)
} # }
```
