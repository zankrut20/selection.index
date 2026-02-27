# Estimate Missing Values in Experimental Data

Estimates and imputes missing values in randomized complete block design
(RCBD), Latin square design (LSD), or split plot design (SPD)
experimental data.

Uses one of six methods: REML, Yates, Healy, Regression, Mean, or
Bartlett.

## Usage

``` r
estimate_missing_values(
  data,
  genotypes,
  replications,
  columns = NULL,
  main_plots = NULL,
  design = c("RCBD", "LSD", "SPD"),
  method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"),
  tolerance = 1e-06
)
```

## Arguments

- data:

  Matrix or data.frame with observations (rows) by traits (columns). May
  contain missing values (NA, NaN, Inf).

- genotypes:

  Vector indicating genotype/treatment for each observation (sub-plot
  treatments in SPD).

- replications:

  Vector indicating replication/block (RCBD) or row (LSD) for each
  observation.

- columns:

  Vector indicating column index for each observation (required for
  Latin Square Design only).

- main_plots:

  Vector indicating main plot treatment for each observation (required
  for Split Plot Design only).

- design:

  Character string specifying experimental design:

  - "RCBD" - Randomized Complete Block Design (default)

  - "LSD" - Latin Square Design

  - "SPD" - Split Plot Design

- method:

  Character string specifying the estimation method:

  - **REML** - Restricted Maximum Likelihood with BLUP. Most robust for
    complex missing patterns. (RCBD, LSD only)

  - **Yates** - Traditional iterative formula. Simple and fast. Good for
    simple missing patterns. (RCBD, LSD only)

  - **Healy** - Healy & Westmacott weighted adjustment method. More
    stable than Yates for multiple missing values. (RCBD, LSD only)

  - **Regression** - Linear regression with QR decomposition.
    Non-iterative, fast and stable. (RCBD, LSD only)

  - **Mean** - Mean substitution using treatment and block effects.
    Non-iterative, fastest. (RCBD, LSD, SPD)

  - **Bartlett** - ANCOVA using other traits as covariates. Best when
    traits are correlated. (RCBD, LSD only)

- tolerance:

  Numeric convergence criterion for iterative methods. Iteration stops
  when maximum change in estimated values falls below this threshold.
  Default: 1e-6.

## Value

Matrix of the same dimensions as `data` with all missing values replaced
by estimates.

## Details

The function handles missing values by iteratively estimating them based
on the experimental design structure:

\*\*RCBD:\*\* 2-way blocking (genotypes × blocks) \*\*LSD:\*\* 3-way
blocking (genotypes × rows × columns) \*\*SPD:\*\* Nested structure
(blocks \> main plots \> sub-plots)

**Method Availability:**

- RCBD: All methods (REML, Yates, Healy, Regression, Mean, Bartlett)

- LSD: All methods (REML, Yates, Healy, Regression, Mean, Bartlett)

- SPD: Mean only (other methods fall back to Mean)

**Method Selection Guide:**

- Use **REML** for complex missing patterns or when precision is
  critical

- Use **Yates** for balanced designs with few missing values

- Use **Healy** when multiple values missing from same treatment/block

- Use **Regression** for fast, deterministic estimation

- Use **Mean** for quick estimation when precision is less critical

- Use **Bartlett** when traits are highly correlated

The function uses the centralized design_stats engine for all ANOVA
computations, ensuring consistency with gen_varcov(), phen_varcov(), and
mean_performance().

## References

Yates, F. (1933). The analysis of replicated experiments when the field
results are incomplete. *Empire Journal of Experimental Agriculture*, 1,
129-142.

Healy, M. J. R., & Westmacott, M. (1956). Missing values in experiments
analysed on automatic computers. *Applied Statistics*, 5(3), 203-206.

Bartlett, M. S. (1937). Some examples of statistical methods of research
in agriculture and applied biology. *Supplement to the Journal of the
Royal Statistical Society*, 4(2), 137-183.

## Examples

``` r
# RCBD example with missing values
data(seldata)
test_data <- seldata[, 3:5]
test_data[c(1, 10, 25), 1] <- NA
test_data[c(5, 15), 2] <- NA

# Impute using Yates method
imputed <- estimate_missing_values(test_data, seldata$treat, seldata$rep, method = "Yates")

# Check that no NA remain
anyNA(imputed) # Should be FALSE
#> [1] FALSE

if (FALSE) { # \dontrun{
# Latin Square Design example
# lsd_data should have genotypes, rows, and columns
imputed_lsd <- estimate_missing_values(
  data = lsd_data[, 3:7],
  genotypes = lsd_data$treat,
  replications = lsd_data$row,
  columns = lsd_data$col,
  design = "LSD",
  method = "REML"
)

# Split Plot Design example
# spd_data should have sub-plots, blocks, and main plots
imputed_spd <- estimate_missing_values(
  data = spd_data[, 3:7],
  genotypes = spd_data$subplot,
  replications = spd_data$block,
  main_plots = spd_data$mainplot,
  design = "SPD",
  method = "Mean"
)
} # }
```
