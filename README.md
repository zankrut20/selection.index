
<!-- README.md is generated from README.Rmd. Please edit that file -->

# selection.index

<!-- badges: start -->

[![CRAN
checks](https://cranchecks.info/badges/summary/selection.index)](https://cran.r-project.org/web/checks/check_results_selection.index.html)
[![R-CMD-check](https://github.com/zankrut20/selection.index/workflows/R-CMD-check/badge.svg)](https://github.com/zankrut20/selection.index/actions)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/zankrut20/selection.index?branch=master&svg=true)](https://ci.appveyor.com/project/zankrut20/selection.index)
[![CRAN
status](https://www.r-pkg.org/badges/version/selection.index)](https://CRAN.R-project.org/package=selection.index)
[![Dependencies](https://tinyverse.netlify.com/badge/selection.index)](https://cran.r-project.org/package=selection.index)
[![Total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/badger?color=blue)](https://cran.r-project.org/package=badger)
[![CRAN
downloads](http://cranlogs.r-pkg.org/badges/last-month/selection.index?color=green)](https://cran.r-project.org/package=selection.index)
[![Weekly
downloads](http://cranlogs.r-pkg.org/badges/last-week/badger?color=yellow)](https://cran.r-project.org/package=badger)
[![coveralls code
coverage](https://coveralls.io/repos/github/google/benchmark/badge.svg?branch=master)](https://coveralls.io/github/google/benchmark)
[![CodeFactor](https://www.codefactor.io/repository/github/zankrut20/selection.index/badge)](https://www.codefactor.io/repository/github/zankrut20/selection.index)
<!-- badges: end -->

The goal of selection.index is to easily construct the selection index
and based on the these indices select the plant traits for the overall
improvement of the plant.

## Installation

You can install the released version of selection.index from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("selection.index")
```

## Example

This is a basic example which shows you how to solve a common problem:
Dataset selindexdata is included in package.

``` r
library(selection.index)
head(seldata)
#>   rep treat   sypp     dtf    rpp    ppr     ppp    spp     pw
#> 1   1    G1 5.4306 42.5000 2.8333 2.0085  7.5833 2.7020 0.5523
#> 2   2    G1 5.4583 42.5000 3.2000 3.7179  7.8000 2.5152 0.7431
#> 3   3    G1 5.5278 43.3333 3.1250 4.2023  7.6111 3.0976 0.7473
#> 4   1    G2 6.3250 43.3333 1.7500 3.0897  3.1000 2.6515 0.4824
#> 5   2    G2 5.8333 43.3333 3.0500 3.7692 14.6500 3.2121 0.6804
#> 6   3    G2 7.9074 43.3333 3.2778 3.6752 12.0000 3.0640 0.6471
```

### Genotypic Variance-Covariance Matrix

``` r
genMat<- gen.varcov(data = seldata[,3:9], genotypes = seldata[,2],
                    replication = seldata[,1])
print(genMat)
#>        sypp     dtf     rpp     ppr     ppp     spp      pw
#> sypp 1.2566  0.3294  0.1588  0.2430  0.7350  0.1276  0.0926
#> dtf  0.3294  1.5602  0.1734 -0.3129 -0.2331  0.1168  0.0330
#> rpp  0.1588  0.1734  0.1325 -0.0316  0.3201 -0.0086 -0.0124
#> ppr  0.2430 -0.3129 -0.0316  0.2432  0.3019 -0.0209  0.0074
#> ppp  0.7350 -0.2331  0.3201  0.3019  0.9608 -0.0692 -0.0582
#> spp  0.1276  0.1168 -0.0086 -0.0209 -0.0692  0.0174  0.0085
#> pw   0.0926  0.0330 -0.0124  0.0074 -0.0582  0.0085  0.0103
```

### Phenotypic Variance-Covariance Matrix

``` r
phenMat<- phen.varcov(data = seldata[,3:9], genotypes = seldata[,2],
                      replication = seldata[,1])
print(phenMat)
#>        sypp     dtf     rpp     ppr     ppp     spp      pw
#> sypp 2.1465  0.1546  0.2320  0.2761  1.0801  0.1460  0.0875
#> dtf  0.1546  3.8372  0.1314 -0.4282 -0.4703  0.0585 -0.0192
#> rpp  0.2320  0.1314  0.2275 -0.0405  0.4635  0.0096 -0.0006
#> ppr  0.2761 -0.4282 -0.0405  0.4678  0.3931 -0.0205  0.0064
#> ppp  1.0801 -0.4703  0.4635  0.3931  4.2638  0.0632 -0.0245
#> spp  0.1460  0.0585  0.0096 -0.0205  0.0632  0.0836  0.0259
#> pw   0.0875 -0.0192 -0.0006  0.0064 -0.0245  0.0259  0.0226
```

### Construction of selection index/indices

For the construction of selection index we requires **phenotypic &
genotypic variance-covariance matrix as well weight matrix.**<br>

-   In **“selection.index”** function **“GAY”** is optional imput. The
    same function is used to calculate the genetic gain of desired
    trait/character’s selection index.

``` r
s<- list()
s[[1]]<- sel.index(ID = 1, phen_mat = phenMat[1,1], 
                   gen_mat = genMat[1,1], weight_mat = weight[1,2])
```

-   Based on above GAY value we calculate the further selection index
    for the other traits/characters.

``` r
s[[2]]<- sel.index(ID = 2, phen_mat = phenMat[2,2],
                   gen_mat = genMat[2,2], weight_mat = weight[2,2], 
                   GAY = 2.1468)
```

### Ranking of selection indices

``` r
r<- rank.index(s,2)
print(r)
#>    ID      b     GA      PRE Rank
#> 2   1 0.5854 1.7694 100.0000    1
#> 21  2 0.4066 1.6431  76.5386    2
```

### Selection score and Ranking of genotypes

``` r
sr<- sel.score.rank(data = seldata[,3], bmat = 0.5, genotype = seldata[,2])
head(sr)
#>   Genotype Selection.score Rank
#> 1       G1        2.736117   19
#> 2       G2        3.344283   13
#> 3       G3        2.276133   23
#> 4       G4        3.503600    9
#> 5       G5        3.506950    8
#> 6       G6        3.068400   17
```
