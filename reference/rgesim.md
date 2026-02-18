# Restricted Linear Genomic Eigen Selection Index Method (RGESIM)

Implements the RGESIM which extends GESIM to allow restrictions on
genetic gains of certain traits. Uses the eigen approach with Lagrange
multipliers.

## Usage

``` r
rgesim(pmat, gmat, Gamma, U_mat, selection_intensity = 2.063)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- Gamma:

  Covariance between phenotypes and GEBVs (n_traits x n_traits).

- U_mat:

  Restriction matrix (r x n_traits) where r is number of restrictions.
  Each row specifies a linear combination of traits to be held at zero
  gain.

- selection_intensity:

  Selection intensity constant \\k_I\\ (default: 2.063 for 10%
  selection).

## Value

Object of class `"rgesim"`, a list with:

- `summary`:

  Data frame with coefficients and metrics.

- `b_y`:

  Coefficients for phenotypic data.

- `b_gamma`:

  Coefficients for GEBVs.

- `b_combined`:

  Combined coefficient vector.

- `E_RG`:

  Expected genetic gains per trait.

- `constrained_response`:

  U' \* E (should be near zero).

- `sigma_I`:

  Standard deviation of the index.

- `hI2`:

  Index heritability.

- `rHI`:

  Accuracy.

- `R_RG`:

  Selection response.

- `lambda2`:

  Leading eigenvalue.

- `implied_w`:

  Implied economic weights.

- `K_RG`:

  Projection matrix.

- `Q_RG`:

  Constraint projection matrix.

- `selection_intensity`:

  Selection intensity used.

## Details

**Eigenproblem (Section 8.4):**
\$\$(\mathbf{K}\_{RG}\mathbf{\Phi}^{-1}\mathbf{A} - \lambda\_{RG}^2
\mathbf{I}\_{2t})\boldsymbol{\beta}\_{RG} = 0\$\$

where: \$\$\mathbf{K}\_{RG} = \mathbf{I}\_{2t} - \mathbf{Q}\_{RG}\$\$
\$\$\mathbf{Q}\_{RG} =
\mathbf{\Phi}^{-1}\mathbf{A}\mathbf{U}\_G(\mathbf{U}\_G^{\prime}\mathbf{A}\mathbf{\Phi}^{-1}\mathbf{A}\mathbf{U}\_G)^{-1}\mathbf{U}\_G^{\prime}\mathbf{A}\$\$

**Implied economic weights:** \$\$\mathbf{w}\_{RG} =
\mathbf{A}^{-1}\[\mathbf{\Phi} +
\mathbf{Q}\_{RG}^{\prime}\mathbf{A}\]\boldsymbol{\beta}\_{RG}\$\$

**Selection response:** \$\$R\_{RG} = k_I
\sqrt{\boldsymbol{\beta}\_{RG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}\_{RG}}\$\$

**Expected genetic gain per trait:** \$\$\mathbf{E}\_{RG} = k_I
\frac{\mathbf{A}\boldsymbol{\beta}\_{RG}}{\sqrt{\boldsymbol{\beta}\_{RG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}\_{RG}}}\$\$

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 8.4.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate GEBV covariance
Gamma <- gmat * 0.8

# Restrict first trait to zero gain
# U_mat must be (n_traits x n_restrictions)
n_traits <- nrow(gmat)
U_mat <- matrix(0, n_traits, 1)
U_mat[1, 1] <- 1  # Restrict trait 1

result <- rgesim(pmat, gmat, Gamma, U_mat)
print(result)
print(result$constrained_response)  # Should be near zero
} # }
```
