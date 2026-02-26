# Getting Started with selection.index

## Introduction

The `selection.index` package provides an optimized suite of tools to
calculate multi-trait phenotypic, genomic, and multistage selection
indices. It simplifies plant and animal breeding workflows by enabling
breeders to accurately rank genotypes based on user-defined economic
weights and experimental designs.

## Setup

First, load the package and the built-in phenotypic dataset. While the
package supports any standard multi-environment phenotypic data, we will
use the included `seldata` dataset (a built-in multi-environment trial
dataset for genotypes) for this quick-start guide.

``` r
library(selection.index)

# Load the built-in phenotypic dataset
data("seldata")

# Inspect the structure of the dataset
head(seldata)
#>   rep treat   sypp     dtf    rpp    ppr     ppp    spp     pw
#> 1   1    G1 5.4306 42.5000 2.8333 2.0085  7.5833 2.7020 0.5523
#> 2   2    G1 5.4583 42.5000 3.2000 3.7179  7.8000 2.5152 0.7431
#> 3   3    G1 5.5278 43.3333 3.1250 4.2023  7.6111 3.0976 0.7473
#> 4   1    G2 6.3250 43.3333 1.7500 3.0897  3.1000 2.6515 0.4824
#> 5   2    G2 5.8333 43.3333 3.0500 3.7692 14.6500 3.2121 0.6804
#> 6   3    G2 7.9074 43.3333 3.2778 3.6752 12.0000 3.0640 0.6471
```

## The Basics

To calculate a selection index, you must specify the traits of interest,
assign economic weights, and compute the genotypic and phenotypic
variance-covariance matrices.

``` r
# Define economic weights for the 7 traits of interest
weights <- c(10, 8, 6, 4, 2, 1, 1)

# Calculate genotypic and phenotypic variance-covariance matrices
# Traits: columns 3:9, Genotypes: column 2, Replication: column 1
gmat <- gen_varcov(data = seldata[, 3:9], genotypes = seldata[, 2], replication = seldata[, 1])
pmat <- phen_varcov(data = seldata[, 3:9], genotypes = seldata[, 2], replication = seldata[, 1])
```

## Execution

With the matrices and weights defined, pass them into the primary
selection function. We use
[`lpsi()`](https://zankrut20.github.io/selection.index/reference/lpsi.md)
to compute the Combinatorial Linear Phenotypic Selection Index for all
traits. Here, setting `ncomb = 7` calculates the index considering all 7
traits simultaneously.

``` r
# Calculate the combinatorial selection index for all 7 traits
index_results <- lpsi(
  ncomb = 7,
  pmat = pmat,
  gmat = gmat,
  wmat = as.matrix(weights),
  wcol = 1
)
```

## Results

Once the index is computed, review the genetic advance metrics and
extract the final selection scores to identify the top-performing
genotypes. `index_results` is a data frame containing the index
coefficients (`b` columns) and various genetic metrics (GA, PRE,
Delta_G, rHI).

``` r
# View the calculated index metrics for our 7-trait combination
head(index_results)
#>                    ID    b.1    b.2    b.3    b.4    b.5    b.6     b.7      GA
#> 1 1, 2, 3, 4, 5, 6, 7 1.0966 1.7309 7.0056 3.2702 0.0596 6.5862 11.5951 21.3233
#>        PRE Delta_G    rHI    hI2 Rank
#> 1 2132.327 21.3233 0.5412 0.2689    1

# Extract the final selection scores to rank the genotypes
scores <- predict_selection_score(
  index_results,
  data = seldata[, 3:9],
  genotypes = seldata[, 2]
)

# View the top ranked genotypes based on their selection scores
head(scores)
#>   Genotypes I_1_2_3_4_5_6_7 I_1_2_3_4_5_6_7_Rank
#> 1        G1        138.8605                   14
#> 2        G2        139.8725                   12
#> 3        G3        134.0505                   21
#> 4        G4        142.2628                    8
#> 5        G5        141.1228                   10
#> 6        G6        136.8680                   18
```
