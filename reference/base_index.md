# Base Index (Williams, 1962)

Implements the Base Index where coefficients are set equal to economic
weights. This is a simple, non-optimized approach that serves as a
baseline comparison.

Unlike the Smith-Hazel index which requires matrix inversion, the Base
Index is computationally trivial and robust when covariance estimates
are unreliable.

## Usage

``` r
base_index(
  pmat,
  gmat,
  wmat,
  wcol = 1,
  selection_intensity = 2.063,
  compare_to_lpsi = TRUE,
  GAY = NULL
)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits)

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits)

- wmat:

  Economic weights matrix (n_traits x k), or vector

- wcol:

  Weight column to use if wmat has multiple columns (default: 1)

- selection_intensity:

  Selection intensity constant (default: 2.063)

- compare_to_lpsi:

  Logical. If TRUE, compares Base Index efficiency to optimal LPSI
  (default: TRUE)

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with:

- `summary` - Data frame with coefficients and metrics

- `b` - Vector of Base Index coefficients (equal to w)

- `w` - Named vector of economic weights

- `Delta_G` - Named vector of expected genetic gains per trait

- `lpsi_comparison` - Optional comparison with Smith-Hazel LPSI

## Details

**Mathematical Formulation:**

Index coefficients: \\b = w\\

The Base Index is appropriate when: - Covariance estimates are
unreliable - Computational simplicity is required - A baseline for
comparison is needed

## Examples

``` r
if (FALSE) { # \dontrun{
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
weights <- c(10, 8, 6, 4, 2, 1, 1)

result <- base_index(pmat, gmat, weights, compare_to_lpsi = TRUE)
print(result)
} # }
```
