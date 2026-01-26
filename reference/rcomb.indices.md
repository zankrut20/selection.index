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
pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
rcomb.indices(ncomb = 2, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,2:3], wcol = 1)
#>      ID    b.1     b.2     GA      PRE Rank
#> 1  2, 3 0.4139  1.1056 2.1197 211.9749    1
#> 2  2, 4 0.3435  0.1655 1.3321 133.2071    7
#> 3  2, 5 0.3718  0.2117 1.6599 165.9936    6
#> 4  2, 6 0.4170  1.3141 1.9305 193.0488    3
#> 5  2, 7 0.4266  2.2775 1.8199 181.9913    5
#> 6  3, 4 0.5321  0.4984 0.8231  82.3055   10
#> 7  3, 5 1.7696  0.1080 1.9995 199.9493    2
#> 8  3, 6 0.5426  0.0426 0.5363  53.6263   13
#> 9  3, 7 0.5279 -0.0766 0.5202  52.0159   14
#> 10 4, 5 0.9935  0.2045 1.8450 184.5050    4
#> 11 4, 6 0.4787  0.0759 0.6722  67.2241   12
#> 12 4, 7 0.5271  0.6326 0.7808  78.0786   11
#> 13 5, 6 0.2208 -0.7864 1.0055 100.5483    9
#> 14 5, 7 0.2007 -1.9031 1.0767 107.6706    8
#> 15 6, 7 0.0807  0.7421 0.2617  26.1727   15
```
