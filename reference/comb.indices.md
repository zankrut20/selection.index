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
pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
wmat<- weight.mat(weight)
comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1, GAY = 1.075)
#>   ID    b.1     GA      PRE Rank
#> 1  1 0.5854 1.7694 164.5979    1
#> 2  2 0.4066 1.6431 152.8479    2
#> 3  3 0.5824 0.5731  53.3070    5
#> 4  4 0.5200 0.7337  68.2467    4
#> 5  5 0.2253 0.9599  89.2920    3
#> 6  6 0.2083 0.1242  11.5579    7
#> 7  7 0.4559 0.1414  13.1535    6
```
