# Phenotypic Variance-Covariance Analysis

Phenotypic Variance-Covariance Analysis

Phenotypic Variance-Covariance Analysis

## Usage

``` r
phen_varcov(
  data,
  genotypes,
  replication,
  columns = NULL,
  main_plots = NULL,
  design_type = c("RCBD", "LSD", "SPD"),
  method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")
)

phen_varcov(
  data,
  genotypes,
  replication,
  columns = NULL,
  main_plots = NULL,
  design_type = c("RCBD", "LSD", "SPD"),
  method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")
)
```

## Arguments

- data:

  traits to be analyzed

- genotypes:

  vector containing genotypes/treatments (sub-plot treatments in SPD)

- replication:

  vector containing replication/blocks (RCBD) or rows (LSD)

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

A Phenotypic Variance-Covariance Matrix

A Phenotypic Variance-Covariance Matrix

## Examples

``` r
# RCBD example
phen_varcov(data=seldata[,3:9], genotypes=seldata$treat, replication=seldata$rep)
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

# Latin Square Design example (requires columns parameter)
# phen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat, 
#            replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")

# Split Plot Design example (requires main_plots parameter)
# phen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot, 
#            replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
# RCBD example
phen_varcov(data=seldata[,3:9], genotypes=seldata$treat, replication=seldata$rep)
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

# Latin Square Design example (requires columns parameter)
# phen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat, 
#            replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")

# Split Plot Design example (requires main_plots parameter)
# phen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot, 
#            replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
```
