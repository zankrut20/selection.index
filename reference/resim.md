# Restricted Linear Phenotypic Eigen Selection Index (RESIM)

Extends ESIM by imposing null restrictions: genetic gains for \\r\\
selected traits are forced to zero while the index heritability for the
remaining traits is maximised.

## Usage

``` r
resim(
  pmat,
  gmat,
  restricted_traits = NULL,
  U_mat = NULL,
  selection_intensity = 2.063
)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- restricted_traits:

  Integer vector of trait indices to restrict to zero genetic gain.
  Example: `c(1, 3)` restricts traits 1 and 3. Alternatively supply a
  custom restriction matrix via `U_mat`.

- U_mat:

  Optional. Restriction matrix (n_traits x r) where each column defines
  one null restriction (\\\mathbf{U}^{\prime}\mathbf{C}\mathbf{b} =
  0\\). Ignored if `restricted_traits` is provided.

- selection_intensity:

  Selection intensity constant (default: 2.063).

## Value

Object of class `"resim"`, a list with:

- `summary`:

  Data frame with b coefficients and key metrics.

- `b`:

  Named numeric vector of RESIM coefficients.

- `Delta_G`:

  Named vector of expected genetic gains per trait.

- `sigma_I`:

  Index standard deviation.

- `hI2`:

  Index heritability (leading eigenvalue of KP^(-1)C).

- `rHI`:

  Index accuracy.

- `lambda2`:

  Leading eigenvalue of the restricted eigenproblem.

- `K`:

  Projection matrix used.

- `U_mat`:

  Restriction matrix used.

- `restricted_traits`:

  Integer vector of restricted trait indices.

- `implied_w`:

  Implied economic weights.

- `selection_intensity`:

  Selection intensity used.

## Details

**Projection matrix (Section 7.2):** \$\$\mathbf{K} = \mathbf{I}\_t -
\mathbf{P}^{-1}\mathbf{C}\mathbf{U}
(\mathbf{U}^{\prime}\mathbf{C}\mathbf{P}^{-1}\mathbf{C}\mathbf{U})^{-1}
\mathbf{U}^{\prime}\mathbf{C}\$\$

**Restricted eigenproblem:** \$\$(\mathbf{K}\mathbf{P}^{-1}\mathbf{C} -
\lambda_R^2 \mathbf{I}\_t)\mathbf{b}\_R = 0\$\$

**Selection response and genetic gain:** \$\$R_R = k_I
\sqrt{\mathbf{b}\_R^{\prime}\mathbf{P}\mathbf{b}\_R}\$\$
\$\$\mathbf{E}\_R = k_I
\frac{\mathbf{C}\mathbf{b}\_R}{\sqrt{\mathbf{b}\_R^{\prime}\mathbf{P}\mathbf{b}\_R}}\$\$

**Implied economic weights:** \$\$\mathbf{w}\_R =
\mathbf{C}^{-1}\[\mathbf{P} +
\mathbf{Q}\_R^{\prime}\mathbf{C}\]\mathbf{b}\_R\$\$ where
\\\mathbf{Q}\_R = \mathbf{I} - \mathbf{K}\\.

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 7.2.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Restrict traits 1 and 3 to zero genetic gain
result <- resim(pmat, gmat, restricted_traits = c(1, 3))
print(result)
summary(result)
} # }
```
