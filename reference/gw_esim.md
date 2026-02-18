# Genome-Wide Linear Eigen Selection Index Method (GW-ESIM)

Implements the GW-ESIM by incorporating genome-wide marker effects
directly into the eigen selection index framework. Uses N marker scores
alongside phenotypic data.

## Usage

``` r
gw_esim(pmat, gmat, G_M, M, selection_intensity = 2.063)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- G_M:

  Covariance between phenotypes and marker scores (n_traits x
  N_markers).

- M:

  Variance-covariance matrix of marker scores (N_markers x N_markers).

- selection_intensity:

  Selection intensity constant \\k_I\\ (default: 2.063 for 10%
  selection).

## Value

Object of class `"gw_esim"`, a list with:

- `summary`:

  Data frame with key metrics.

- `b_y`:

  Coefficients for phenotypic data.

- `b_m`:

  Coefficients for marker scores.

- `b_combined`:

  Combined coefficient vector.

- `E_W`:

  Expected genetic gains per trait.

- `sigma_I`:

  Standard deviation of the index.

- `hI2`:

  Index heritability (= leading eigenvalue).

- `rHI`:

  Accuracy.

- `R_W`:

  Selection response.

- `lambda2`:

  Leading eigenvalue.

- `selection_intensity`:

  Selection intensity used.

## Details

**Eigenproblem (Section 8.3):** \$\$(\mathbf{Q}^{-1}\mathbf{X} -
\lambda_W^2 \mathbf{I}\_{(t+N)})\boldsymbol{\beta}\_W = 0\$\$

where: \$\$\mathbf{Q} = \begin{bmatrix} \mathbf{P} &
\mathbf{G}\_M^{\prime} \\ \mathbf{G}\_M & \mathbf{M} \end{bmatrix}\$\$
\$\$\mathbf{X} = \begin{bmatrix} \mathbf{C} & \mathbf{G}\_M^{\prime} \\
\mathbf{G}\_M & \mathbf{M} \end{bmatrix}\$\$

**Selection response:** \$\$R_W = k_I
\sqrt{\boldsymbol{\beta}\_W^{\prime}\mathbf{Q}\boldsymbol{\beta}\_W}\$\$

**Expected genetic gain per trait:** \$\$\mathbf{E}\_W = k_I
\frac{\mathbf{X}\boldsymbol{\beta}\_W}{\sqrt{\boldsymbol{\beta}\_W^{\prime}\mathbf{Q}\boldsymbol{\beta}\_W}}\$\$

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 8.3.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate marker data
N_markers <- 100
n_traits <- nrow(gmat)
G_M <- matrix(rnorm(n_traits * N_markers, sd = 0.5), n_traits, N_markers)
M <- diag(N_markers) + matrix(rnorm(N_markers^2, sd = 0.1), N_markers, N_markers)
M <- (M + t(M)) / 2  # Make symmetric

result <- gw_esim(pmat, gmat, G_M, M)
print(result)
} # }
```
