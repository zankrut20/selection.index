# Genotypic Variance-Covariance Analysis

Genotypic Variance-Covariance Analysis

## Usage

``` r
gen_varcov(
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

A Genotypic Variance-Covariance Matrix

## Examples

``` r
# RCBD example
gen_varcov(data = seldata[, 3:9], genotypes = seldata$treat, replication = seldata$rep)
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

# Latin Square Design example (requires columns parameter)
# gen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat,
#           replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")

# Split Plot Design example (requires main_plots parameter)
# gen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot,
#           replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
```
