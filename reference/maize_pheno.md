# Synthetic Maize Phenotypic Data

A synthetic dataset containing multi-environment phenotypic records for
100 maize genotypes. Designed to demonstrate phenotypic selection index
calculations and variance-covariance matrix estimations.

## Usage

``` r
data(maize_pheno)
```

## Format

A data frame with 600 rows and 6 variables:

- Genotype:

  Factor representing 100 unique maize genotypes.

- Environment:

  Factor representing 2 distinct growing environments.

- Block:

  Factor representing 3 replicate blocks within each environment.

- Yield:

  Numeric vector of grain yield in kg/ha.

- PlantHeight:

  Numeric vector of plant height in centimeters.

- DaysToMaturity:

  Numeric vector of days to physiological maturity.

## Source

Simulated for the selection.index package to provide reproducible
examples.

## Examples

``` r
data(maize_pheno)
head(maize_pheno)
#>   Genotype Environment Block   Yield PlantHeight DaysToMaturity
#> 1     G001        Env1    B1 8260.29      224.83            117
#> 2     G002        Env1    B1 7460.15      197.42            119
#> 3     G003        Env1    B1 8069.62      211.43            117
#> 4     G004        Env1    B1 7957.63      214.27            116
#> 5     G005        Env1    B1 7666.68      198.77            119
#> 6     G006        Env1    B1 6741.96      182.73            116
```
