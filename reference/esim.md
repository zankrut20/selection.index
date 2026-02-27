# Linear Phenotypic Eigen Selection Index (ESIM)

Implements the ESIM by maximising the squared accuracy \\\rho\_{HI}^2\\
through the generalised eigenproblem of the multi-trait heritability
matrix \\\mathbf{P}^{-1}\mathbf{C}\\.

Unlike the Smith-Hazel LPSI, \*\*no economic weights are required\*\*.
The net genetic merit vector \\\mathbf{w}\_E\\ is instead implied by the
solution.

## Usage

``` r
esim(pmat, gmat, selection_intensity = 2.063, n_indices = 1L)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).
  Corresponds to **C** in the Chapter 7 notation.

- selection_intensity:

  Selection intensity constant \\k_I\\ (default: 2.063 for 10%
  selection).

- n_indices:

  Number of leading ESIM vectors to return (default: 1). Returning \>1
  provides a ranked set of indices for comparative analysis.

## Value

Object of class `"esim"`, a list with:

- `summary`:

  Data frame with b coefficients, hI2, rHI, sigma_I, Delta_G, and
  lambda2 for each index requested.

- `b`:

  Named numeric vector of optimal ESIM coefficients (1st index).

- `Delta_G`:

  Named numeric vector of expected genetic gains per trait.

- `sigma_I`:

  Standard deviation of the index \\\sigma_I\\.

- `hI2`:

  Index heritability \\h^2\_{I_E}\\ (= leading eigenvalue).

- `rHI`:

  Accuracy \\r\_{HI_E}\\.

- `lambda2`:

  Leading eigenvalue (maximised index heritability).

- `implied_w`:

  Implied economic weights \\\mathbf{w}\_E\\.

- `all_eigenvalues`:

  All eigenvalues of \\\mathbf{P}^{-1}\mathbf{C}\\.

- `selection_intensity`:

  Selection intensity used.

## Details

**Eigenproblem (Section 7.1):** \$\$(\mathbf{P}^{-1}\mathbf{C} -
\lambda_E^2 \mathbf{I})\mathbf{b}\_E = 0\$\$

The solution \\\lambda_E^2\\ (largest eigenvalue) equals the maximum
achievable index heritability \\h^2\_{I_E}\\.

**Key metrics:** \$\$R_E = k_I
\sqrt{\mathbf{b}\_E^{\prime}\mathbf{P}\mathbf{b}\_E}\$\$
\$\$\mathbf{E}\_E = k_I
\frac{\mathbf{C}\mathbf{b}\_E}{\sqrt{\mathbf{b}\_E^{\prime}\mathbf{P}\mathbf{b}\_E}}\$\$

**Implied economic weights:** \$\$\mathbf{w}\_E =
\frac{\sqrt{\lambda_E^2}}{\mathbf{b}\_E^{\prime}\mathbf{P}\mathbf{b}\_E}
\mathbf{C}^{-1}\mathbf{P}\mathbf{b}\_E\$\$

Uses `cpp_symmetric_solve` and `cpp_quadratic_form_sym` from
`math_primitives.cpp` for efficient matrix operations, and R's
[`eigen()`](https://rdrr.io/r/base/eigen.html) for the
eigendecomposition.

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 7.1.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

result <- esim(pmat, gmat)
print(result)
summary(result)
} # }
```
