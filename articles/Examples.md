# Data Analysis with selection.index

The aim of most plant breeding program is simultaneous improvement of
several characters. An objective method involving simultaneous selection
for several attributes then becomes necessary. It has been recognized
that most rapid improvements in the economic value is expected from
selection applied simultaneously to all the characters which determine
the economic value of a plant, and appropriate assigned weights to each
character according to their economic importance, heritability and
correlations between characters. So the selection for economic value is
a complex matter. If the component characters are combined together into
an index in such a way that when selection is applied to the index, as
if index is the character to be improved, most rapid improvement of
economic value is excepted. Such an index was first proposed by Smith
(1937) based on the Fisher’s (1936) “discriminant function”. In this
package selection index is calculated based on the Smith (1937)
selection index method (Dabholkar, 1999). For more information refer
**Elements of Bio Metrical GENETICS by A. R. Dabholkar.**

``` r
library(selection.index)
d<- seldata # Manually generated data for analysis which is included in package
```

``` r
w<- weight # Weights assigned to the traits also include in package
```

As we discussed that selection index based on discriminant function. So
we have required **genotypic & phenotypic variance-covariance matrix**
for further analysis.

- Genotypic variance-covariance matrix

``` r
gmat<- gen_varcov(data = d[,3:9], genotypes = d$treat, replication = d$rep)
print(gmat)
#>            sypp         dtf          rpp          ppr         ppp          spp
#> sypp 1.25660210  0.32936305  0.158785900  0.242981986  0.73499020  0.127571993
#> dtf  0.32936305  1.56017847  0.173388420 -0.312908175 -0.23310004  0.116790239
#> rpp  0.15878590  0.17338842  0.132484364 -0.031596521  0.32014873 -0.008643769
#> ppr  0.24298199 -0.31290818 -0.031596521  0.243231727  0.30192365 -0.020860985
#> ppp  0.73499020 -0.23310004  0.320148725  0.301923650  0.96076644 -0.069172364
#> spp  0.12757199  0.11679024 -0.008643769 -0.020860985 -0.06917236  0.017410958
#> pw   0.09261588  0.03298807 -0.012353519  0.007352443 -0.05824420  0.008560105
#>                pw
#> sypp  0.092615879
#> dtf   0.032988075
#> rpp  -0.012353519
#> ppr   0.007352443
#> ppp  -0.058244197
#> spp   0.008560105
#> pw    0.010304709
```

- Phenotypic variance-covariance matrix

``` r
pmat<- phen_varcov(data = d[,3:9], genotypes = d$treat, replication = d$rep)
print(pmat)
#>           sypp         dtf          rpp         ppr        ppp          spp
#> sypp 4.6596932  0.81327830  0.549544688  0.76208582  2.5500613  0.401129893
#> dtf  0.8132783  6.95753031  0.478114232 -1.05398088 -0.9364611  0.292048297
#> rpp  0.5495447  0.47811423  0.492447900 -0.10364377  1.1038273 -0.007695489
#> ppr  0.7620858 -1.05398088 -0.103643767  0.95426114  0.9969895 -0.062186773
#> ppp  2.5500613 -0.93646109  1.103827327  0.99698955  6.1852816 -0.075103699
#> spp  0.4011299  0.29204830 -0.007695489 -0.06218677 -0.0751037  0.118394770
#> pw   0.2727074  0.04678408 -0.025308347  0.02107724 -0.1410159  0.043030640
#>               pw
#> sypp  0.27270744
#> dtf   0.04678408
#> rpp  -0.02530835
#> ppr   0.02107724
#> ppp  -0.14101586
#> spp   0.04303064
#> pw    0.04321272
```

Generally, **Percent Relative Efficiency (PRE)** of a selection index is
calculated with reference to **Genetic Advance (GA) yield** of
respective weight. So first we calculate the GA of yield for respective
weights. + Genetic gain of Yield

``` r
GAY<- gen_advance(phen_mat = pmat[1,1], gen_mat = gmat[1,1],
                  weight_mat = w[1,2])
print(GAY)
#>        [,1]
#> [1,] 1.2009
```

We use this GAY value for the construction, ranking of the other
selection indices and stored them in a list “si”.

## Selection score and Ranking of genotypes

Generally selection score is calculate based on top ranked selection
index. So first we store the **discriminant coefficient** value into a
variable **b**, and later that value we used for calculation of
selection score and ranking of the genotypes.

## `lpsi()` is used for construction of selection indices based on different combination of characters.

``` r
lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)
#>   ID    b.1     GA      PRE Delta_G    rHI    hI2 Rank
#> 1  1 0.6316 2.8125 234.2021  2.8125 0.4825 0.2697    1
#> 2  2 0.2396 1.3036 108.5481  1.3036 0.2237 0.2242    4
#> 3  3 1.4869 2.1526 179.2462  2.1526 0.3693 0.2690    2
#> 4  4 0.4507 0.9084  75.6402  0.9084 0.1558 0.2549    6
#> 5  5 0.3164 1.6236 135.1988  1.6236 0.2786 0.1553    3
#> 6  6 1.4499 1.0292  85.7009  1.0292 0.1766 0.1471    5
#> 7  7 1.8796 0.8061  67.1225  0.8061 0.1383 0.2385    7
```

## Using `excluding_trait` parameter to remove traits from the construction of selection indices

``` r
lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY, excluding_trait = 1)
#>   ID    b.1     GA      PRE Delta_G    rHI    hI2 Rank
#> 1  2 0.2396 1.3036 108.5481  1.3036 0.2237 0.2242    3
#> 2  3 1.4869 2.1526 179.2462  2.1526 0.3693 0.2690    1
#> 3  4 0.4507 0.9084  75.6402  0.9084 0.1558 0.2549    5
#> 4  5 0.3164 1.6236 135.1988  1.6236 0.2786 0.1553    2
#> 5  6 1.4499 1.0292  85.7009  1.0292 0.1766 0.1471    4
#> 6  7 1.8796 0.8061  67.1225  0.8061 0.1383 0.2385    6
```
