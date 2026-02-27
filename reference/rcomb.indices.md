# Remove trait or trait combination from possible trait combinations of possible Trait combinations

Remove trait or trait combination from possible trait combinations of
possible Trait combinations

## Usage

``` r
rcomb.indices(ncomb, i, pmat, gmat, wmat, wcol = 1, GAY)
```

## Arguments

- ncomb:

  Number of character combination

- i:

  remove trait or trait combination

- pmat:

  Phenotypic Variance Covariance Matrix

- gmat:

  Genotypic Variance Covariance Matrix

- wmat:

  Weight Matrix

- wcol:

  Respective weight column number of Weight Matrix

- GAY:

  Genetic Advance/Genetic Gain of base selection index

## Value

Data frame of possible selection indices with per cent relative
efficiency and ranking

## Examples

``` r
gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#> Error in gen.varcov(seldata[, 3:9], seldata[, 2], seldata[, 1]): could not find function "gen.varcov"
pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#> Error in phen.varcov(seldata[, 3:9], seldata[, 2], seldata[, 1]): could not find function "phen.varcov"
rcomb.indices(ncomb = 2, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,2:3], wcol = 1)
#> Error in rcomb.indices(ncomb = 2, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,     2:3], wcol = 1): could not find function "rcomb.indices"
```
