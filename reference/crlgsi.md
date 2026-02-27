# Combined Restricted Linear Genomic Selection Index (CRLGSI)

Implements the CRLGSI which combines phenotypic and genomic information
with restrictions on genetic gains. This extends CLGSI to include
constraints.

## Usage

``` r
crlgsi(
  T_C = NULL,
  Psi_C = NULL,
  phen_mat = NULL,
  gebv_mat = NULL,
  pmat = NULL,
  gmat = NULL,
  wmat,
  wcol = 1,
  restricted_traits = NULL,
  U = NULL,
  reliability = NULL,
  k_I = 2.063,
  L_I = 1,
  GAY = NULL
)
```

## Arguments

- T_C:

  Combined variance-covariance matrix (2t x 2t) where t = n_traits.
  Structure: \[P, P_yg; P_yg', P_g\] where P = phenotypic var, P_g =
  GEBV var, P_yg = covariance between phenotypes and GEBVs. Can be
  computed automatically if phen_mat and gebv_mat are provided.

- Psi_C:

  Combined genetic covariance matrix (2t x t). Structure: \[G;
  C_gebv_g\] where G = genetic var, C_gebv_g = Cov(GEBV, g). Can be
  computed automatically if gmat and reliability are provided.

- phen_mat:

  Optional. Matrix of phenotypes (n_genotypes x n_traits)

- gebv_mat:

  Optional. Matrix of GEBVs (n_genotypes x n_traits)

- pmat:

  Optional. Phenotypic variance-covariance matrix

- gmat:

  Optional. Genotypic variance-covariance matrix

- wmat:

  Economic weights matrix (n_traits x k), or vector

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- restricted_traits:

  Vector of trait indices to restrict (default: NULL)

- U:

  Constraint matrix (2t x n_constraints for combined traits).
  Alternative to restricted_traits. Ignored if restricted_traits is
  provided.

- reliability:

  Optional. Reliability of GEBVs (r^2)

- k_I:

  Selection intensity (default: 2.063)

- L_I:

  Standardization constant (default: 1)

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with:

- `summary` - Data frame with coefficients and metrics

- `b` - Vector of CRLGSI coefficients (β_CR)

- `b_y` - Coefficients for phenotypes

- `b_g` - Coefficients for GEBVs

- `E` - Expected genetic gains per trait

- `R` - Overall selection response

## Details

**Mathematical Formulation (Chapter 6, Section 6.3):**

The CRLGSI combines phenotypic and genomic data with restrictions.

Coefficient vector: \\beta_CR = K_C \* beta_C\\

Where K_C incorporates the restriction matrix.

Selection response: \\R_CR = (k_I / L_I) \* sqrt(beta_CR' \* T_C \*
beta_CR)\\

Expected gains: \\E_CR = (k_I / L_I) \* (Psi_C \* beta_CR) /
sqrt(beta_CR' \* T_C \* beta_CR)\\

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
n_genotypes <- 100
n_traits <- 5

phen_mat <- matrix(rnorm(n_genotypes * n_traits, 15, 3), n_genotypes, n_traits)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, 10, 2), n_genotypes, n_traits)

gmat <- cov(phen_mat) * 0.6 # Genotypic component
pmat <- cov(phen_mat)

w <- c(10, 8, 6, 4, 2)

# Restrict traits 2 and 4
result <- crlgsi(
  phen_mat = phen_mat, gebv_mat = gebv_mat,
  pmat = pmat, gmat = gmat, wmat = w,
  restricted_traits = c(2, 4), reliability = 0.7
)
print(result$summary)
} # }
```
