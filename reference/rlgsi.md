# Restricted Linear Genomic Selection Index (RLGSI)

Implements the Restricted Linear Genomic Selection Index where genetic
gains are constrained to zero for specific traits while maximizing gains
for others. Uses GEBVs only (no phenotypic data required).

## Usage

``` r
rlgsi(
  Gamma,
  wmat,
  wcol = 1,
  restricted_traits = NULL,
  U = NULL,
  k_I = 2.063,
  L_G = 1,
  gmat = NULL,
  GAY = NULL
)
```

## Arguments

- Gamma:

  GEBV variance-covariance matrix (n_traits x n_traits). This represents
  the variance of GEBVs, typically computed from predicted breeding
  values.

- wmat:

  Economic weights matrix (n_traits x k), or vector

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- restricted_traits:

  Vector of trait indices to restrict (default: NULL). Example: c(1, 3)
  restricts traits 1 and 3 to zero gain.

- U:

  Constraint matrix (n_traits x n_constraints). Each column defines a
  restriction. Alternative to restricted_traits for custom constraints.
  Ignored if restricted_traits is provided.

- k_I:

  Selection intensity (default: 2.063 for 10 percent selection)

- L_G:

  Standardization constant (default: 1). Can be set to sqrt(w'Gw) for
  standardization.

- gmat:

  Optional. True genetic variance-covariance matrix for exact accuracy
  calculation. If NULL, uses Gamma as approximation. Providing gmat
  ensures textbook-perfect accuracy metric.

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with:

- `summary` - Data frame with coefficients, response metrics

- `b` - Vector of RLGSI coefficients (\\\beta\\\_RG)

- `E` - Named vector of expected genetic gains per trait

- `R` - Overall selection response

- `U` - Constraint matrix used

- `constrained_response` - Realized gains for constrained traits (should
  be ~0)

## Details

**Mathematical Formulation (Chapter 6, Section 6.1):**

The RLGSI minimizes the mean squared difference between the index I =
\\\beta\\'\\\gamma\\ and the breeding objective H = w'g under the
restriction: U'\\\Gamma\\\\\beta\\ = 0

Solution involves solving the augmented system: \$\$\begin{bmatrix}
\Gamma & \Gamma U \\ U'\Gamma & 0 \end{bmatrix} \begin{bmatrix} \beta \\
v \end{bmatrix} = \begin{bmatrix} \Gamma w \\ 0 \end{bmatrix}\$\$

Where: - \\\Gamma\\ (Gamma) = Var(GEBVs) - GEBV variance-covariance
matrix - U = Constraint matrix (each column is a restriction vector) - w
= Economic weights - \\\beta\\\_RG = RLGSI coefficient vector - v =
Lagrange multipliers

Selection response: \\R\_{RG} = (k_I / L_G) \* sqrt(beta_RG' \* Gamma \*
beta_RG)\\

Expected gains: \\E\_{RG} = (k_I / L_G) \* (Gamma \* beta_RG) /
sqrt(beta_RG' \* Gamma \* beta_RG)\\

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate GEBV variance-covariance matrix
set.seed(123)
n_traits <- 5
Gamma <- matrix(rnorm(n_traits^2), n_traits, n_traits)
Gamma <- (Gamma + t(Gamma)) / 2 # Make symmetric
diag(Gamma) <- abs(diag(Gamma)) + 2 # Ensure positive definite

# Economic weights
w <- c(10, 8, 6, 4, 2)

# Restrict traits 2 and 4 to zero gain
result <- rlgsi(Gamma, w, restricted_traits = c(2, 4))
print(result$summary)
print(result$E) # Check that traits 2 and 4 have ~0 gain
} # }
```
