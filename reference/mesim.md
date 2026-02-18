# Molecular Eigen Selection Index Method (MESIM)

Implements the MESIM by maximising the squared accuracy through the
generalised eigenproblem combining phenotypic data with molecular marker
scores. Unlike Smith-Hazel LPSI, \*\*no economic weights are
required\*\*.

## Usage

``` r
mesim(pmat, gmat, S_M, S_Mg = NULL, S_var = NULL, selection_intensity = 2.063)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- S_M:

  Covariance between phenotypes and marker scores (n_traits x n_traits).
  This represents Cov(y, s) where s are marker scores. Used in the
  phenotypic variance matrix T_M.

- S_Mg:

  Covariance between genetic values and marker scores (n_traits x
  n_traits). This represents Cov(g, s). Used in the genetic variance
  matrix Psi_M. If NULL (default), uses S_M, which is appropriate when
  assuming Cov(e, s) ~= 0 (errors uncorrelated with markers), so
  Cov(y,s) ~= Cov(g,s).

- S_var:

  Variance-covariance matrix of marker scores (n_traits x n_traits).
  Represents Var(s). If NULL (default), uses S_M for backward
  compatibility, but providing the actual Var(s) is more statistically
  rigorous.

- selection_intensity:

  Selection intensity constant \\k_I\\ (default: 2.063 for 10%
  selection).

## Value

Object of class `"mesim"`, a list with:

- `summary`:

  Data frame with coefficients and metrics.

- `b_y`:

  Coefficients for phenotypic data.

- `b_s`:

  Coefficients for marker scores.

- `b_combined`:

  Combined coefficient vector.

- `E_M`:

  Expected genetic gains per trait.

- `sigma_I`:

  Standard deviation of the index.

- `hI2`:

  Index heritability (= leading eigenvalue).

- `rHI`:

  Accuracy \\r\_{HI}\\.

- `R_M`:

  Selection response.

- `lambda2`:

  Leading eigenvalue (maximised index heritability).

- `selection_intensity`:

  Selection intensity used.

## Details

**Eigenproblem (Section 8.1):**
\$\$(\mathbf{T}\_M^{-1}\mathbf{\Psi}\_M - \lambda_M^2
\mathbf{I}\_{2t})\boldsymbol{\beta}\_M = 0\$\$

where: \$\$\mathbf{T}\_M = \begin{bmatrix} \mathbf{P} &
\mathrm{Cov}(\mathbf{y},\mathbf{s}) \\
\mathrm{Cov}(\mathbf{s},\mathbf{y}) & \mathrm{Var}(\mathbf{s})
\end{bmatrix}\$\$ \$\$\mathbf{\Psi}\_M = \begin{bmatrix} \mathbf{C} &
\mathrm{Cov}(\mathbf{g},\mathbf{s}) \\
\mathrm{Cov}(\mathbf{s},\mathbf{g}) & \mathrm{Var}(\mathbf{s})
\end{bmatrix}\$\$

**Theoretical distinction:**

- T_M uses phenotypic covariances: Cov(y, s) provided via `S_M`

- Psi_M uses genetic covariances: Cov(g, s) provided via `S_Mg`

- Since y = g + e, if Cov(e, s) ~= 0, then Cov(y, s) ~= Cov(g, s)

- Chapter 8.1 assumes marker scores are pure genetic predictors, so S_M
  can be used for both (default behavior when S_Mg = NULL)

The solution \\\lambda_M^2\\ (largest eigenvalue) equals the maximum
achievable index heritability.

**Key metrics:** \$\$R_M = k_I
\sqrt{\boldsymbol{\beta}\_{M}^{\prime}\mathbf{T}\_M\boldsymbol{\beta}\_{M}}\$\$
\$\$\mathbf{E}\_M = k_I
\frac{\mathbf{\Psi}\_M\boldsymbol{\beta}\_{M}}{\sqrt{\boldsymbol{\beta}\_{M}^{\prime}\mathbf{T}\_M\boldsymbol{\beta}\_{M}}}\$\$

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 8.1.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate marker score matrices (in practice, compute from data)
S_M <- gmat * 0.7      # Cov(y, s) - phenotype-marker covariance
S_Mg <- gmat * 0.65    # Cov(g, s) - genetic-marker covariance
S_var <- gmat * 0.8    # Var(s) - marker score variance

# Most rigorous: Provide all three covariance matrices
result <- mesim(pmat, gmat, S_M, S_Mg = S_Mg, S_var = S_var)
print(result)

# Standard usage: Cov(g,s) defaults to Cov(y,s) when errors uncorrelated
result_standard <- mesim(pmat, gmat, S_M, S_var = S_var)

# Backward compatible: Chapter 8.1 simplified notation
result_simple <- mesim(pmat, gmat, S_M)
} # }
```
