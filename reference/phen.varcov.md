# Phenotypic Variance-Covariance Analysis

Phenotypic Variance-Covariance Analysis

## Usage

``` r
phen.varcov(data, genotypes, replication)
```

## Arguments

- data:

  traits to be analyzed

- genotypes:

  vector containing genotypes/treatments

- replication:

  vector containing replication

## Value

A Phenotypic Variance-Covariance Matrix

## Examples

``` r
phen.varcov(data=seldata[,3:9], genotypes=seldata$treat,replication=seldata$rep)
#>            sypp         dtf           rpp          ppr         ppp          spp
#> sypp 2.14648906  0.15455221  0.2319728887  0.276121850  1.08008088  0.145985907
#> dtf  0.15455221  3.83717336  0.1313373906 -0.428164534 -0.47026101  0.058467819
#> rpp  0.23197289  0.13133739  0.2274791728 -0.040450725  0.46352988  0.009592048
#> ppr  0.27612185 -0.42816453 -0.0404507247  0.467797686  0.39314225 -0.020464804
#> ppp  1.08008088 -0.47026101  0.4635298765  0.393142248  4.26374874  0.063241030
#> spp  0.14598591  0.05846782  0.0095920480 -0.020464804  0.06324103  0.083572855
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
