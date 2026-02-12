# Construction of selection indices (filtered/combinatorial)

Build all possible selection indices from trait combinations, with
optional exclusion of specific traits (by index, name, or data columns).

## Usage

``` r
lpsi(ncomb, pmat, gmat, wmat, wcol = 1, GAY, excluding_trait = NULL)
```

## Arguments

- ncomb:

  Number of Characters/Traits group

- pmat:

  Phenotypic Variance-Covariance Matrix

- gmat:

  Genotypic Variance-Covariance Matrix

- wmat:

  Weight Matrix

- wcol:

  Weight column number incase more than one weights, by default its 1

- GAY:

  Genetic Advance of comparative Character/Trait i.e. Yield (Optional
  argument)

- excluding_trait:

  Optional. Traits to exclude from combinations. Can be: (1) numeric
  vector of trait indices (e.g., c(1, 3)), (2) character vector of trait
  names (e.g., c("sypp", "dtf")), (3) data frame/matrix columns with
  trait data (trait names extracted from column names). When specified,
  only combinations that do NOT contain any of these traits are
  returned.

## Value

Data frame of all possible selection indices with metrics (GA, PRE,
Delta_G, rHI, hI2)

## Examples

``` r
gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
wmat<- weight_mat(weight)
lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1, GAY = 1.075)
#>   ID    b.1     GA      PRE Delta_G    rHI    hI2 Rank
#> 1  1 0.2697 1.2009 111.7146  1.2009 0.5193 0.2697    2
#> 2  2 0.2242 1.2202 113.5109  1.2202 0.4735 0.2242    1
#> 3  3 0.2690 0.3895  36.2306  0.3895 0.5187 0.2690    5
#> 4  4 0.2549 0.5137  47.7834  0.5137 0.5049 0.2549    4
#> 5  5 0.1553 0.7970  74.1359  0.7970 0.3941 0.1553    3
#> 6  6 0.1471 0.1044   9.7106  0.1044 0.3835 0.1471    6
#> 7  7 0.2385 0.1023   9.5131  0.1023 0.4883 0.2385    7
```
