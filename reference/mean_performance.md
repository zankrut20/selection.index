# Mean performance of phenotypic data

Mean performance of phenotypic data

## Usage

``` r
mean_performance(
  data,
  genotypes,
  replications,
  columns = NULL,
  main_plots = NULL,
  design_type = c("RCBD", "LSD", "SPD"),
  method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")
)
```

## Arguments

- data:

  data for analysis

- genotypes:

  genotypes vector (sub-plot treatments in SPD)

- replications:

  replication vector

- columns:

  vector containing columns (required for Latin Square Design only)

- main_plots:

  vector containing main plot treatments (required for Split Plot Design
  only)

- design_type:

  experimental design type: "RCBD" (default), "LSD" (Latin Square), or
  "SPD" (Split Plot)

- method:

  Method for missing value imputation: "REML" (default), "Yates",
  "Healy", "Regression", "Mean", or "Bartlett"

## Value

Dataframe of mean performance analysis

## Examples

``` r
mean_performance(data = seldata[, 3:9], genotypes = seldata[, 2], replications = seldata[, 1])
#>          Genotypes    sypp     dtf     rpp     ppr       ppp       spp      pw
#> 1               G1  5.4722 42.7778  3.0528  3.3096    7.6648    2.7716  0.6809
#> 2               G2  6.6886 43.3333  2.6926  3.5114    9.9167    2.9759  0.6033
#> 3               G3  4.5523 45.0000  1.9405  2.4756    5.3135     3.189  0.7036
#> 4               G4  7.0072 44.1667  2.5500  3.3149    9.1233    3.2364  0.6526
#> 5               G5  7.0139 43.8889  2.5517  3.4179    9.0833    3.0849  0.6513
#> 6               G6  6.1368 43.0555  2.0667  2.8366    6.6375    3.3055  0.8359
#> 7               G7  4.9012 43.0555  2.1806  2.6361    5.9431    2.9613  0.6761
#> 8               G8  5.5278 45.8333  2.4417  3.1168    7.6875     2.883  0.5705
#> 9               G9  8.5787 43.0556  2.2639  3.9815    7.4282    3.0415  0.9469
#> 10             G10  3.8491 43.8889  1.9444  3.0413    5.8222     2.633  0.6218
#> 11             G11  7.0748 44.7222  2.1620  3.4081    6.4491    3.1748  0.9267
#> 12             G12  6.6061 43.0555  1.7379  3.7300    5.9402    3.0557  0.8817
#> 13             G13  6.7768 43.3333  1.9926  3.5000    6.3241    3.0185  0.8398
#> 14             G14  8.1632 43.0555  2.6083  4.2528    8.4208    3.0795  0.7400
#> 15             G15  4.2917 44.4444  1.3056  3.9732    5.4815    2.9307  0.6044
#> 16             G16  6.5927 43.6111  2.7672  3.4045    9.4709    3.0198  0.6834
#> 17             G17  6.9761 43.6111  2.5000  3.1718    8.9793     2.997  0.6490
#> 18             G18  6.2523 43.0555  2.4924  4.0041    8.5348    2.9679  0.6435
#> 19             G19  7.7929 46.6667  2.5204  2.5698    6.6426    3.5118  0.9076
#> 20             G20  7.6460 42.2222  1.6481  4.2279    7.0185     3.459  0.9427
#> 21             G21  7.8182 45.8333  2.3352  3.1918    8.1148    3.2166  0.7387
#> 22             G22  7.6370 49.1666  2.8333  2.9459    6.3555     3.266  0.8327
#> 23             G23  5.2694 43.8889  2.3185  2.8056    5.9037    3.2192  0.6252
#> 24             G24  5.4680 43.0555  1.9778  2.3801    4.8652    3.0723  0.6928
#> 25             G25  6.7222 45.2778  2.5741  2.3148    6.4852    3.2884  0.8429
#> 26             Min  3.8491 42.2222  1.3056  2.3148    4.8652     2.633  0.5705
#> 27             Max  8.5787 49.1666  3.0528  4.2528    9.9167    3.5118  0.9469
#> 28              GM  6.4326 44.1222  2.2983  3.2609    7.1843    3.0944  0.7398
#> 29          CV (%) 14.6650  3.4200 13.4103 14.5324   25.2971    8.3125 14.9911
#> 30             SEm  0.5446  0.8712  0.1779  0.2736    1.0493    0.1485  0.0640
#> 31           CD 5%  1.5487  2.4772  0.5060  0.7780    2.9836    0.4223  0.1821
#> 32           CD 1%  2.0659  3.3047  0.6750  1.0378 3.9801 NS 0.5633 NS  0.2429
#> 33    Heritability  0.5854  0.4066  0.5824  0.5200    0.2253    0.2083  0.4559
#> 34 Heritability(%) 58.5422 40.6596 58.2402 51.9951   22.5334   20.8333 45.5894
```
