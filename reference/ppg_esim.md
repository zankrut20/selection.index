# Predetermined Proportional Gain Eigen Selection Index (PPG-ESIM)

Extends ESIM by enforcing that genetic gains are proportional to a
user-specified vector \\\mathbf{d}\\: \\\Delta\mathbf{G} \propto
\mathbf{d}\\. A similarity transformation \\\boldsymbol{\beta}\_P =
\mathbf{F}\mathbf{b}\_P\\ aligns the eigenvector with the breeder's
desired direction.

## Usage

``` r
ppg_esim(pmat, gmat, d, selection_intensity = 2.063)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits).

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits).

- d:

  Numeric vector of desired proportional gains (length n_traits). The
  *ratios* among elements define target gain proportions. Direction
  (positive/negative) must reflect desired improvement direction
  (positive = increase, negative = decrease).

- selection_intensity:

  Selection intensity constant (default: 2.063).

## Value

Object of class `"ppg_esim"`, a list with:

- `summary`:

  Data frame with beta (transformed b), b (raw), hI2, rHI, sigma_I,
  Delta_G, and lambda2.

- `beta`:

  Named numeric vector of post-transformation PPG-ESIM coefficients
  \\\boldsymbol{\beta}\_P = \mathbf{F}\mathbf{b}\_P\\.

- `b`:

  Raw eigenvector b_P before similarity transform.

- `Delta_G`:

  Named vector of expected genetic gains per trait.

- `sigma_I`:

  Index standard deviation.

- `hI2`:

  Index heritability.

- `rHI`:

  Index accuracy.

- `lambda2`:

  Leading eigenvalue of the PPG restricted eigenproblem.

- `F_mat`:

  Diagonal similarity transform matrix F (diag(sign(d))).

- `K_P`:

  PPG projection matrix (rank 1: projects onto d subspace).

- `D_M`:

  Mallard matrix (t x t-1): orthogonal complement of d, used to
  construct the (t-1) restrictions.

- `desired_gains`:

  Input proportional gains vector d.

- `selection_intensity`:

  Selection intensity used.

## Details

**Restriction structure via the Mallard Matrix (Section 7.3):**

The PPG-ESIM restricts the \\(t-1)\\ directions \*\*orthogonal\*\* to
\\\mathbf{d}\\, forcing the genetic gain vector to be collinear with
\\\mathbf{d}\\.

The Mallard matrix \\\mathbf{D}\_M\\ is \\t \times (t-1)\\: its columns
span the orthogonal complement of \\\mathbf{d}\\, obtained via QR
decomposition of \\\mathbf{d}/\\\mathbf{d}\\\\: \$\$\mathbf{Q}\_{QR} =
\[\hat{d} \mid \mathbf{D}\_M\], \quad \text{QR}(\hat{d}) \to
\mathbf{Q}\_{QR} \in \mathbb{R}^{t \times t}\$\$

With \\\boldsymbol{\Psi} = \mathbf{C}\\ (full-trait case, \\\mathbf{U} =
\mathbf{I}\_t\\):

**PPG projection matrix (\\t-1\\ restrictions):** \$\$\mathbf{Q}\_P =
\mathbf{P}^{-1}\boldsymbol{\Psi}\mathbf{D}\_M
(\mathbf{D}\_M^{\prime}\boldsymbol{\Psi}^{\prime}\mathbf{P}^{-1}\boldsymbol{\Psi}\mathbf{D}\_M)^{-1}
\mathbf{D}\_M^{\prime}\boldsymbol{\Psi}^{\prime}\$\$ \$\$\mathbf{K}\_P =
\mathbf{I}\_t - \mathbf{Q}\_P \quad (\text{rank 1})\$\$

Because \\\mathbf{K}\_P\\ has rank 1 (projects onto the \\\mathbf{d}\\
subspace), \\\mathbf{K}\_P\mathbf{P}^{-1}\mathbf{C}\\ has exactly one
positive eigenvalue and its eigenvector produces \\\Delta\mathbf{G}
\propto \mathbf{d}\\.

**PPG eigenproblem (rank-1 system):**
\$\$(\mathbf{K}\_P\mathbf{P}^{-1}\mathbf{C} -
\lambda_P^2\mathbf{I}\_t)\mathbf{b}\_P = 0\$\$

**Similarity transform:** \$\$\boldsymbol{\beta}\_P =
\mathbf{F}\mathbf{b}\_P\$\$ where \\\mathbf{F} =
\text{diag}(\text{sign}(\mathbf{d}))\\ aligns the eigenvector sign with
the breeder's intended improvement direction.

**Key response metrics:** \$\$R_P =
k_I\sqrt{\boldsymbol{\beta}\_P^{\prime}\mathbf{P}\boldsymbol{\beta}\_P}\$\$
\$\$\mathbf{E}\_P =
k_I\frac{\mathbf{C}\boldsymbol{\beta}\_P}{\sqrt{\boldsymbol{\beta}\_P^{\prime}\mathbf{P}\boldsymbol{\beta}\_P}}\$\$

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Section 7.3.

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Desired proportional gains: increase all traits proportionally
d <- c(2, 1, 1, 1, 1, 1, 1)
result <- ppg_esim(pmat, gmat, d)
print(result)
summary(result)
} # }
```
