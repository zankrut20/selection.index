# The Linear Phenotypic Selection Index Theory

## Introduction

In plant and animal breeding, quantitative traits (QTs) are expressions
of genes distributed across the genome interacting with the environment.
The phenotypic value of QTs ($y$) can be systematically partitioned into
a genotypic component ($g$) and an environmental component ($e$):

$$y = g + e$$

The primary goal in breeding is to maximize an individual’s **net
genetic merit**. The net genetic merit ($H$) is a linear combination of
the unobservable true breeding values ($\mathbf{g}$) weighted by their
respective economic values ($\mathbf{w}$):

$$H = \mathbf{w}^{\prime}\mathbf{g}$$

Because the net genetic merit is unobservable in field trials, breeders
construct a **Linear Phenotypic Selection Index (LPSI)** to predict it.
The LPSI ($I$) is a linear combination of the observable and optimally
weighted phenotypic trait values ($\mathbf{y}$) adjusted by index
coefficients ($\mathbf{b}$):

$$I = \mathbf{b}^{\prime}\mathbf{y}$$

The objective of the LPSI is to predict the net genetic merit and
maximize the multi-trait selection response.

## Optimizing the LPSI

To identify the optimal parents for the next selection cycle, the
correlation between the net genetic merit ($H$) and the LPSI ($I$) must
be maximized. The vector $\mathbf{b}$ that simultaneously minimizes the
mean squared difference between $I$ and $H$ and perfectly maximizes this
correlation is mathematically derived as:

$$\mathbf{b} = \mathbf{P}^{- 1}{\mathbf{G}\mathbf{w}}$$

where: \* $\mathbf{P}$ is the phenotypic variance-covariance matrix. \*
$\mathbf{G}$ is the genotypic variance-covariance matrix. \*
$\mathbf{w}$ is the vector of economic weights defining relative trait
importance.

Once these optimal coefficients are derived, we can evaluate two
fundamental parameters:

1.  **The Maximized Selection Response ($R_{I}$)**: The expected mean
    improvement in the net genetic merit due to indirect selection on
    the index.
    $$R_{I} = k_{I}\sqrt{\mathbf{b}^{\prime}{\mathbf{P}\mathbf{b}}}$$

2.  **The Expected Genetic Gain Per Trait ($\mathbf{E}$)**: The
    multi-trait selection response broken down per individual trait.
    $$\mathbf{E} = k_{I}\frac{\mathbf{G}\mathbf{b}}{\sigma_{I}}$$

where $k_{I}$ is the standardized selection intensity and $\sigma_{I}$
is the standard deviation of the index score variance.

## Practical Implementation in R

We can seamlessly translate this text theory into rigorous statistical
practice using the `selection.index` package. We will utilize the
built-in synthetic datasets: `maize_pheno` (containing multi-environment
phenotypic records for 100 genotypes) and `maize_geno` (500 SNP
markers).

### 1. Estimating Covariance Matrices

First, we estimate the genotypic ($\mathbf{G}$) and phenotypic
($\mathbf{P}$) variance-covariance matrices from our raw phenotypic
dataset.

``` r
library(selection.index)

# Load the synthetic phenotypic multi-environment dataset
data("maize_pheno")

# In maize_pheno: Traits are columns 4:6.
# Genotypes are in column 1, and Block/Replication is in column 3.
gmat <- gen_varcov(data = maize_pheno[, 4:6], genotypes = maize_pheno[, 1], replication = maize_pheno[, 3])
pmat <- phen_varcov(data = maize_pheno[, 4:6], genotypes = maize_pheno[, 1], replication = maize_pheno[, 3])
```

### 2. Defining Economic Weights

Next, we establish the relative economic priority of each trait.
Economic weights ($\mathbf{w}$) explicitly define our strategic breeding
objectives.

``` r
# Define the economic weights for the 3 continuous traits
# (e.g., Yield, PlantHeight, DaysToMaturity)
weights <- c(10, -5, -5)
```

### 3. Calculating the LPSI

With the covariance matrices and economic weights specified, we
integrate them into the primary
[`lpsi()`](https://zankrut20.github.io/selection.index/reference/lpsi.md)
function, which evaluates the combinatorial multi-trait selection
indices efficiently.

``` r
# Calculate the Optimal Combinatorial Linear Phenotypic Selection Index (LPSI)
index_results <- lpsi(
  ncomb = 3,
  pmat = pmat,
  gmat = gmat,
  wmat = as.matrix(weights),
  wcol = 1
)
```

### 4. Evaluating Outcomes and Selecting Genotypes

Finally, we evaluate the theoretical gains. The
[`lpsi()`](https://zankrut20.github.io/selection.index/reference/lpsi.md)
function returns a structured data frame containing the theoretical
selection response ($R_{I}$) and other parameter estimates for all
requested trait combinations.

``` r
# View the top combinatorial indices, including their selection response (R_A)
head(index_results)
#>        ID      b.1      b.2     b.3      GA       PRE Delta_G     rHI      hI2
#> 1 1, 2, 3 2364.671 9453.933 1012156 6863519 686351898 6863519 42.0899 1771.794
#>   Rank
#> 1    1

# Extract the phenotypic selection scores to strategically rank the parental candidates
# using the top evaluated combinatorial index
scores <- predict_selection_score(
  index_results,
  data = maize_pheno[, 4:6],
  genotypes = maize_pheno[, 1]
)

# View the top performing candidates designated for the next breeding cycle
head(scores)
#>   Genotypes   I_1_2_3 I_1_2_3_Rank
#> 1      G001 138639943           78
#> 2      G002 138393284           83
#> 3      G003 138260546           88
#> 4      G004 139442740           46
#> 5      G005 139319528           51
#> 6      G006 139661606           36
```

### 5. Extension: Linear Marker Selection Index

The classical linear selection index theories seamlessly extend to
marker-assisted genomic selection. If you have genome-wide marker
profiles for your genotypes, you can incorporate them to estimate the
Linear Marker Selection Index (LMSI).

``` r
# Load the associated synthetic genomic dataset (500 SNPs for the 100 genotypes)
data("maize_geno")

# Calculate the marker-assisted index combining our matrices and raw SNP profiles
marker_index_results <- lmsi(
  pmat = pmat,
  gmat = gmat,
  marker_scores = maize_geno,
  wmat = weights
)

summary(marker_index_results)
```

### 6. The Base Index and Index Efficiency

In scenarios where the phenotypic ($\mathbf{P}$) and genotypic
($\mathbf{G}$) matrices are poorly estimated (e.g., due to limited
data), the true optimal coefficients ($\mathbf{b}$) can be
systematically biased. The **Base Index** provides a robust,
non-optimized alternative where coefficients are set strictly equal to
the fixed economic weights ($I_{B} = \mathbf{w}\prime\mathbf{y}$).

``` r
# Calculate the Base Index and automatically compare its efficiency to the LPSI
base_results <- base_index(
  pmat = pmat,
  gmat = gmat,
  wmat = weights,
  compare_to_lpsi = TRUE
)

# Observe the expected genetic gains and efficiency comparison
base_results$summary
#>   b.1 b.2 b.3      GA       PRE  Delta_G      hI2     rHI
#> 1  10  -5  -5 1824121 182412094 14577.53 125.1324 11.1863
```

### 7. Heritability of the LPSI

The theory demonstrates that the correlation between the net genetic
merit ($H$) and the expected index ($I$) differs from the traditional
index heritability mathematically ($h_{I}^{2} \neq \rho_{HI}^{2}$). The
[`lpsi()`](https://zankrut20.github.io/selection.index/reference/lpsi.md)
function intrinsically estimates both of these fundamental statistics:

``` r
# Extract the top combinatorial index results
top_index <- index_results[1, ]

# h^2_I: Heritability of the optimal index
top_index$hI2
#> [1] 1771.794

# \rho_HI: Correlation between the LPSI and the true underlying Net Genetic Merit
top_index$rHI
#> [1] 42.0899
```
