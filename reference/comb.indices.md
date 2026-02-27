# Construction of selection indices based on number of character grouping

Construction of selection indices based on number of character grouping

## Usage

``` r
comb.indices(ncomb, pmat, gmat, wmat, wcol = 1, GAY)
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

## Value

Data frame of all possible selection indices

## Examples

``` r
gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#> Error in gen.varcov(seldata[, 3:9], seldata[, 2], seldata[, 1]): could not find function "gen.varcov"
pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#> Error in phen.varcov(seldata[, 3:9], seldata[, 2], seldata[, 1]): could not find function "phen.varcov"
wmat<- weight.mat(weight)
#> Error in weight.mat(weight): could not find function "weight.mat"
comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1, GAY = 1.075)
#> Error in comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = wmat,     wcol = 1, GAY = 1.075): could not find function "comb.indices"
```
