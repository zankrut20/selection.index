# Smith-Hazel Linear Phenotypic Selection Index

Implements the optimal Smith-Hazel selection index which maximizes the
correlation between the index I = b'y and the breeding objective H =
w'g.

This is the foundational selection index method from Chapter 2.

## Usage

``` r
smith_hazel(
  pmat,
  gmat,
  wmat,
  wcol = 1,
  selection_intensity = 2.063,
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

  Selection intensity constant (default: 2.063 for 10% selection)

- GAY:

  Optional. Genetic advance of comparative trait for PRE calculation

## Value

List with:

- `summary` - Data frame with coefficients and metrics

- `b` - Vector of Smith-Hazel index coefficients

- `w` - Named vector of economic weights

- `Delta_G` - Named vector of expected genetic gains per trait

- `sigma_I` - Standard deviation of the index

- `GA` - Total genetic advance

- `PRE` - Percent relative efficiency

- `hI2` - Heritability of the index

- `rHI` - Accuracy (correlation with breeding objective)

## Details

**Mathematical Formulation (Chapter 2):**

Index coefficients: \\b = P^{-1}Gw\\

Where: - P = Phenotypic variance-covariance matrix - G = Genotypic
variance-covariance matrix - w = Economic weights

Key metrics: - Variance of index: \\\sigma^2_I = b'Pb\\ - Total genetic
advance: \\R_H = i\sqrt{b'Pb}\\ - Expected gains per trait: \\\Delta G =
(i/\sigma_I)Gb\\ - Heritability of index: \\h^2_I = b'Gb / b'Pb\\ -
Accuracy: \\r\_{HI} = \sqrt{b'Gb / b'Pb}\\

## Examples

``` r
if (FALSE) { # \dontrun{
# Calculate variance-covariance matrices
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Define economic weights
weights <- c(10, 8, 6, 4, 2, 1, 1)

# Build Smith-Hazel index
result <- smith_hazel(pmat, gmat, weights)
print(result)
summary(result)
} # }
```
