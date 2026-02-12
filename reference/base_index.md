# Base Index (Williams, 1962)

Implements the Base Index where selection coefficients equal economic
weights (b = w). This is a simplified index used when reliable estimates
of variance-covariance matrices are unavailable, or as a baseline for
comparison with optimized indices. It assumes zero correlations between
traits but still calculates expected response using the true genetic
covariances.

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

  Phenotypic variance-covariance matrix (p x p)

- gmat:

  Genotypic variance-covariance matrix (p x p)

- wmat:

  Matrix of economic weights (p x k), where k is number of weight sets.
  Can also be a numeric vector which will be converted to a matrix.

- wcol:

  Column index of wmat to use if wmat has multiple columns (default: 1)

- selection_intensity:

  Selection intensity (default: 2.063, corresponding to 10% selection)

- compare_to_lpsi:

  Logical; if TRUE, also calculate and compare to standard LPSI
  (default: TRUE)

- GAY:

  Genetic advance of comparative trait (optional, for PRE calculation)

## Value

List of class base_index with components:

- `b` - Numeric vector of index coefficients (equal to economic weights)

- `w` - Economic weights used

- `Delta_G` - Named vector of expected genetic response per trait

- `sigma_I` - Standard deviation of the index

- `GA` - Genetic advance in the index

- `PRE` - Percent relative efficiency (if GAY provided)

- `hI2` - Heritability of the index

- `rHI` - Correlation between index and aggregate genotype

- `selection_intensity` - Selection intensity used

- `summary` - Data frame with coefficients and metrics

- `lpsi_comparison` - Comparison with LPSI (if compare_to_lpsi = TRUE)

## Details

The Base Index (Williams, 1962) is the simplest selection index where:
\$\$\mathbf{b} = \mathbf{w}\$\$

The expected response is: \$\$\Delta\mathbf{G}\_{base} =
\frac{i}{\sigma_I} \mathbf{G} \mathbf{w}\$\$ where \$\$\sigma_I =
\sqrt{\mathbf{w}' \mathbf{P} \mathbf{w}}\$\$

This index does not optimize the selection response like the Smith-Hazel
LPSI, but provides a simple baseline and is robust when covariance
estimates are unreliable.

## References

Williams, J.S. (1962). The evaluation of a selection index. Biometrics,
18, 375-393.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
data(seldata)

# Calculate variance-covariance matrices
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Define economic weights (e.g., favor first trait)
weights <- c(10, 5, 3, 3, 5, 8, 4)

# Calculate Base Index
result <- base_index(pmat, gmat, weights)
print(result)

# Compare efficiency with LPSI
summary(result)
} # }
```
