# Mean performance of phenotypic data

Mean performance of phenotypic data

## Usage

``` r
meanPerformance(data, genotypes, replications)
```

## Arguments

- data:

  data for analysis

- genotypes:

  genotypes vector

- replications:

  replication vector

## Value

Dataframe of mean performance analysis

## Examples

``` r
meanPerformance(data = seldata[, 3:9], genotypes = seldata[, 2], replications = seldata[, 1])
#> Error in meanPerformance(data = seldata[, 3:9], genotypes = seldata[,     2], replications = seldata[, 1]): could not find function "meanPerformance"

```
