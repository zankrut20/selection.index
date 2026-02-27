# Genotypic Variance-Covariance Analysis

Genotypic Variance-Covariance Analysis

## Usage

``` r
gen.varcov(data, genotypes, replication)
```

## Arguments

- data:

  traits to be analyzed

- genotypes:

  vector containing genotypes/treatments

- replication:

  vector containing replication

## Value

A Genotypic Variance-Covariance Matrix

## Examples

``` r
gen.varcov(data=seldata[,3:9], genotypes=seldata$treat,replication=seldata$rep)
#> Error in gen.varcov(data = seldata[, 3:9], genotypes = seldata$treat,     replication = seldata$rep): could not find function "gen.varcov"
```
