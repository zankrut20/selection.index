# Predetermined Proportional Gains (PPG-LPSI)

Predetermined Proportional Gains (PPG-LPSI)

## Usage

``` r
ppg_lpsi(pmat, gmat, k, wmat = NULL, wcol = 1, GAY)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix

- gmat:

  Genotypic variance-covariance matrix

- k:

  Vector of desired proportional gains

- wmat:

  Optional weight matrix for GA/PRE calculation

- wcol:

  Weight column number (default: 1)

- GAY:

  Genetic advance of comparative trait (optional)

## Value

List with summary data frame, coefficient vector, Delta_G vector, and
phi

## Examples

``` r
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
k <- rep(1, ncol(pmat))
ppg_lpsi(pmat, gmat, k)
#> $summary
#>       b.1    b.2     b.3      b.4    b.5       b.6     b.7 GA PRE Delta_G rHI
#> 1 11.4784 3.8899 -17.521 -17.3052 0.5096 -115.2781 71.4102 NA  NA 53.5871  NA
#>       hI2    phi
#> 1 -0.0931 0.0794
#> 
#> $b
#> [1]   11.4784    3.8899  -17.5210  -17.3052    0.5096 -115.2781   71.4102
#> 
#> $Delta_G
#>       sypp        dtf        rpp        ppr        ppp        spp         pw 
#> 0.07942154 0.07942154 0.07942154 0.07942154 0.07942154 0.07942154 0.07942154 
#> 
#> $phi
#> [1] 0.07942154
#> 
```
