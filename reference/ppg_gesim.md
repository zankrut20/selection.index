# Predetermined Proportional Gain Genomic Eigen Selection Index (PPG-GESIM)

Implements the PPG-GESIM which extends GESIM to enforce that genetic
gains are proportional to a user-specified vector d. Combines eigen
approach with predetermined gain proportions.

## Usage

``` r
ppg_gesim(pmat, gmat, Gamma, d, selection_intensity = 2.063)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- Gamma:

  Covariance between phenotypes and GEBVs (n_traits x n_traits).

- d:

  Numeric vector of desired proportional gains (length n_traits). The
  ratios among elements define target gain proportions.

- selection_intensity:

  Selection intensity constant \\k_I\\ (default: 2.063 for 10%
  selection).

## Value

Object of class `"ppg_gesim"`, a list with:

- `summary`:

  Data frame with coefficients and metrics.

- `b_y`:

  Coefficients for phenotypic data.

- `b_gamma`:

  Coefficients for GEBVs.

- `b_combined`:

  Combined coefficient vector.

- `E_PG`:

  Expected genetic gains per trait.

- `gain_ratios`:

  Ratios of actual to desired gains (should be constant).

- `d`:

  Original desired proportional gains (length t).

- `d_PG`:

  Extended proportional gains (length 2t, includes GEBV targets).

- `sigma_I`:

  Standard deviation of the index.

- `hI2`:

  Index heritability.

- `rHI`:

  Accuracy.

- `R_PG`:

  Selection response.

- `lambda2`:

  Leading eigenvalue.

- `implied_w`:

  Implied economic weights.

- `U_PG`:

  Restriction matrix ((2t-1) x 2t).

- `selection_intensity`:

  Selection intensity used.

## Details

**Eigenproblem (Section 8.5):** \$\$(\mathbf{T}\_{PG} - \lambda\_{PG}^2
\mathbf{I}\_{2t})\boldsymbol{\beta}\_{PG} = 0\$\$

where: \$\$\mathbf{T}\_{PG} =
\mathbf{K}\_{RG}\mathbf{\Phi}^{-1}\mathbf{A} + \mathbf{B}\$\$
\$\$\mathbf{B} = \boldsymbol{\delta}\boldsymbol{\varphi}^{\prime}\$\$

**Implied economic weights:** \$\$\mathbf{w}\_{PG} =
\mathbf{A}^{-1}\[\mathbf{\Phi} +
\mathbf{Q}\_{PG}^{\prime}\mathbf{A}\]\boldsymbol{\beta}\_{PG}\$\$

**Selection response:** \$\$R\_{PG} = k_I
\sqrt{\boldsymbol{\beta}\_{PG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}\_{PG}}\$\$

**Expected genetic gain per trait:** \$\$\mathbf{E}\_{PG} = k_I
\frac{\mathbf{A}\boldsymbol{\beta}\_{PG}}{\sqrt{\boldsymbol{\beta}\_{PG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}\_{PG}}}\$\$

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 8.5.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Simulate GEBV covariance
Gamma <- gmat * 0.8

# Desired proportional gains (e.g., 2:1:3 ratio for first 3 traits)
d <- c(2, 1, 3, 1, 1, 1, 1)

result <- ppg_gesim(pmat, gmat, Gamma, d)
print(result)
print(result$gain_ratios) # Should be approximately constant
} # }
```
