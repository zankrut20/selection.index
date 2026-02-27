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
#> Error in phen.varcov(data = seldata[, 3:9], genotypes = seldata$treat,     replication = seldata$rep): could not find function "phen.varcov"
```
