
<!-- README.md is generated from README.Rmd. Please edit that file -->

# selection.index

<!-- badges: start -->

[![R-CMD-check](https://github.com/zankrut20/selection.index/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zankrut20/selection.index/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/selection.index)](https://CRAN.R-project.org/package=selection.index)
[![CRAN
checks](https://badges.cranchecks.info/summary/selection.index.svg)](https://cran.r-project.org/web/checks/check_results_selection.index.html)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/selection.index?color=blue)](https://cran.r-project.org/package=selection.index)
[![Monthly
downloads](https://cranlogs.r-pkg.org/badges/last-month/selection.index?color=green)](https://cran.r-project.org/package=selection.index)
[![Codecov test
coverage](https://codecov.io/gh/zankrut20/selection.index/branch/master/graph/badge.svg)](https://app.codecov.io/gh/zankrut20/selection.index?branch=master)
[![CodeFactor](https://www.codefactor.io/repository/github/zankrut20/selection.index/badge)](https://www.codefactor.io/repository/github/zankrut20/selection.index)
[![GitHub
stars](https://img.shields.io/github/stars/zankrut20/selection.index?style=social)](https://github.com/zankrut20/selection.index/stargazers)
[![GitHub
issues](https://img.shields.io/github/issues/zankrut20/selection.index)](https://github.com/zankrut20/selection.index/issues)
[![GitHub
forks](https://img.shields.io/github/forks/zankrut20/selection.index?style=social)](https://github.com/zankrut20/selection.index/network/members)
[![Last
commit](https://img.shields.io/github/last-commit/zankrut20/selection.index)](https://github.com/zankrut20/selection.index/commits/master)
[![R
version](https://img.shields.io/badge/R-%E2%89%A5%203.5.0-blue)](https://www.r-project.org/)

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

from [github](https://github.com/zankrut20/selection.index) with:

``` r
devtools::install_github("zankrut20/selection.index")
```

## Example

This is a basic example which shows you how to solve a common problem:
Dataset `seldata` is included in package.

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
#>            sypp         dtf          rpp          ppr         ppp          spp
#> sypp 1.25660210  0.32936305  0.158785900  0.242981986  0.73499020  0.127571993
#> dtf  0.32936305  1.56017847  0.173388420 -0.312908175 -0.23310004  0.116790239
#> rpp  0.15878590  0.17338842  0.132484364 -0.031596521  0.32014873 -0.008643769
#> ppr  0.24298199 -0.31290818 -0.031596521  0.243231727  0.30192365 -0.020860985
#> ppp  0.73499020 -0.23310004  0.320148725  0.301923650  0.96076644 -0.069172364
#> spp  0.12757199  0.11679024 -0.008643769 -0.020860985 -0.06917236  0.017410958
#> pw   0.09261588  0.03298807 -0.012353519  0.007352443 -0.05824420  0.008560105
#>                pw
#> sypp  0.092615879
#> dtf   0.032988075
#> rpp  -0.012353519
#> ppr   0.007352443
#> ppp  -0.058244197
#> spp   0.008560105
#> pw    0.010304709
```

### Phenotypic Variance-Covariance Matrix

``` r
phenMat<- phen.varcov(data = seldata[,3:9], genotypes = seldata[,2],
                      replication = seldata[,1])
print(phenMat)
#>            sypp         dtf           rpp          ppr         ppp          spp
#> sypp 2.14648906  0.15455221  0.2319728887  0.276121850  1.08008088  0.145985907
#> dtf  0.15455221  3.83717336  0.1313373906 -0.428164534 -0.47026101  0.058467819
#> rpp  0.23197289  0.13133739  0.2274791728 -0.040450725  0.46352988  0.009592048
#> ppr  0.27612185 -0.42816453 -0.0404507247  0.467797686  0.39314225 -0.020464804
#> ppp  1.08008088 -0.47026101  0.4635298765  0.393142248  4.26374874  0.063241030
#> spp  0.14598591  0.05846782  0.0095920481 -0.020464804  0.06324103  0.083572855
#> pw   0.08747568 -0.01919207 -0.0006013091  0.006372357 -0.02452747  0.025910429
#>                 pw
#> sypp  0.0874756848
#> dtf  -0.0191920675
#> rpp  -0.0006013091
#> ppr   0.0063723573
#> ppp  -0.0245274713
#> spp   0.0259104291
#> pw    0.0226032987
```

### Weight Matrix - Data is included in package `weight`

``` r
weightMat <- weight.mat(weight)
weightMat
#>      ew     h2
#> [1,]  1 0.6947
#> [2,]  1 0.3244
#> [3,]  1 0.6861
#> [4,]  1 0.7097
#> [5,]  1 0.8336
#> [6,]  1 0.2201
#> [7,]  1 0.5226
```

- Genetic gain of Yield

``` r
GAY<- gen.advance(phen_mat = phenMat[1,1], gen_mat = genMat[1,1],
                  weight_mat = weightMat[1,1])
print(GAY)
#>        [,1]
#> [1,] 1.7694
```

### Construction of selection index/indices

For the construction of selection index we requires **phenotypic &
genotypic variance-covariance matrix as well weight matrix.**<br>

### Construction of all possible selection indices for a character combinations

``` r
comb.indices(ncomb = 1, pmat = phenMat, gmat = genMat, wmat = weight[,2:3], wcol = 1, GAY = GAY)
#>   ID    b.1     GA      PRE Rank
#> 1  1 0.5854 1.7694 100.0015    1
#> 2  2 0.4066 1.6431  92.8628    2
#> 3  3 0.5824 0.5731  32.3867    5
#> 4  4 0.5200 0.7337  41.4634    4
#> 5  5 0.2253 0.9599  54.2494    3
#> 6  6 0.2083 0.1242   7.0220    7
#> 7  7 0.4559 0.1414   7.9914    6
```

### Construction of selection indices by removing desired single character from the combinations

``` r
rcomb.indices(ncomb = 1, i = 1, pmat = phenMat, gmat = genMat, wmat = weight[,2:3], wcol = 1, GAY = GAY)
#>   ID    b.1     GA     PRE Rank
#> 1  2 0.4066 1.6431 92.8628    1
#> 2  3 0.5824 0.5731 32.3867    4
#> 3  4 0.5200 0.7337 41.4634    3
#> 4  5 0.2253 0.9599 54.2494    2
#> 5  6 0.2083 0.1242  7.0220    6
#> 6  7 0.4559 0.1414  7.9914    5
```
