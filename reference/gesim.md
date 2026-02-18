# Linear Genomic Eigen Selection Index Method (GESIM)

Implements the GESIM by maximising the squared accuracy through the
generalised eigenproblem combining phenotypic data with GEBVs (Genomic
Estimated Breeding Values). No economic weights are required.

## Usage

``` r
gesim(pmat, gmat, Gamma, selection_intensity = 2.063)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- Gamma:

  Covariance between phenotypes and GEBVs (n_traits x n_traits). This
  represents Cov(y, gamma) where gamma are GEBVs.

- selection_intensity:

  Selection intensity constant \\k_I\\ (default: 2.063 for 10%
  selection).

## Value

Object of class `"gesim"`, a list with:

- `summary`:

  Data frame with coefficients and metrics.

- `b_y`:

  Coefficients for phenotypic data.

- `b_gamma`:

  Coefficients for GEBVs.

- `b_combined`:

  Combined coefficient vector.

- `E_G`:

  Expected genetic gains per trait.

- `sigma_I`:

  Standard deviation of the index.

- `hI2`:

  Index heritability (= leading eigenvalue).

- `rHI`:

  Accuracy \\r\_{HI}\\.

- `R_G`:

  Selection response.

- `lambda2`:

  Leading eigenvalue.

- `implied_w`:

  Implied economic weights.

- `selection_intensity`:

  Selection intensity used.

## Details

**Eigenproblem (Section 8.2):** \$\$(\mathbf{\Phi}^{-1}\mathbf{A} -
\lambda_G^2 \mathbf{I}\_{2t})\boldsymbol{\beta}\_G = 0\$\$

where: \$\$\mathbf{\Phi} = \begin{bmatrix} \mathbf{P} & \mathbf{\Gamma}
\\ \mathbf{\Gamma} & \mathbf{\Gamma} \end{bmatrix}\$\$ \$\$\mathbf{A} =
\begin{bmatrix} \mathbf{C} & \mathbf{\Gamma} \\ \mathbf{\Gamma} &
\mathbf{\Gamma} \end{bmatrix}\$\$

**Implied economic weights:** \$\$\mathbf{w}\_G =
\mathbf{A}^{-1}\mathbf{\Phi}\boldsymbol{\beta}\$\$

**Selection response:** \$\$R_G = k_I
\sqrt{\boldsymbol{\beta}\_G^{\prime}\mathbf{\Phi}\boldsymbol{\beta}\_G}\$\$

**Expected genetic gain per trait:** \$\$\mathbf{E}\_G = k_I
\frac{\mathbf{A}\boldsymbol{\beta}\_G}{\sqrt{\boldsymbol{\beta}\_G^{\prime}\mathbf{\Phi}\boldsymbol{\beta}\_G}}\$\$

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 8.2.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate GEBV covariance (in practice, compute from genomic predictions)
Gamma <- gmat * 0.8  # Assume 80% GEBV-phenotype covariance

result <- gesim(pmat, gmat, Gamma)
print(result)
} # }
```
