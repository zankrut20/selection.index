# Genetic Advance for PRE

Genetic Advance for PRE

## Usage

``` r
gen.advance(phen_mat, gen_mat, weight_mat)
```

## Arguments

- phen_mat:

  phenotypic matrix value of desired characters

- gen_mat:

  genotypic matrix value of desired characters

- weight_mat:

  weight matrix value of desired characters

## Value

Genetic advance of character or character combinations

## Examples

``` r
gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
gen.advance(phen_mat = pmat[1,1], gen_mat = gmat[1,1], weight_mat = weight[1,2])
#>        [,1]
#> [1,] 1.7694
```
