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
gmat<- gen.varcov(data = d[,3:9], genotypes = d$treat, replication = d$rep)
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
pmat<- phen.varcov(data = d[,3:9], genotypes = d$treat, replication = d$rep)
print(pmat)
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

Generally, **Percent Relative Efficiency (PRE)** of a selection index is
calculated with reference to **Genetic Advance (GA) yield** of
respective weight. So first we calculate the GA of yield for respective
weights. + Genetic gain of Yield

``` r
GAY<- gen.advance(phen_mat = pmat[1,1], gen_mat = gmat[1,1],
                  weight_mat = w[1,2])
print(GAY)
#>        [,1]
#> [1,] 1.7694
```

We use this GAY value for the construction, ranking of the other
selection indices and stored them in a list “si”.

## Selection score and Ranking of genotypes

Generally selection score is calculate based on top ranked selection
index. So first we store the **discriminant coefficient** value into a
variable **b**, and later that value we used for calculation of
selection score and ranking of the genotypes.

## `comb.indices()` is used for construction of selection indices based on different combination of characters.

``` r
comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)
#>   ID    b.1     GA      PRE Rank
#> 1  1 0.5854 1.7694 100.0015    1
#> 2  2 0.4066 1.6431  92.8628    2
#> 3  3 0.5824 0.5731  32.3867    5
#> 4  4 0.5200 0.7337  41.4634    4
#> 5  5 0.2253 0.9599  54.2494    3
#> 6  6 0.2083 0.1242   7.0220    7
#> 7  7 0.4559 0.1414   7.9914    6
```

## \`rcomb.indices()\`\` - remove trait from the construction of selection indices

``` r
rcomb.indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = w[,-1], wcol = 1, GAY = GAY)
#>   ID    b.1     GA     PRE Rank
#> 1  2 0.4066 1.6431 92.8628    1
#> 2  3 0.5824 0.5731 32.3867    4
#> 3  4 0.5200 0.7337 41.4634    3
#> 4  5 0.2253 0.9599 54.2494    2
#> 5  6 0.2083 0.1242  7.0220    6
#> 6  7 0.4559 0.1414  7.9914    5
```
