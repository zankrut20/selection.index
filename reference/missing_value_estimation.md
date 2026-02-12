# Missing Value Estimation for RCBD, Latin Square, and Split Plot Designs

Estimates missing values in Randomized Complete Block Design (RCBD),
Latin Square Design (LSD), or Split Plot Design (SPD) using one of six
methods: REML, Yates, Healy, Regression, Mean, or Bartlett.

## Usage

``` r
missing_value_estimation(
  data_mat,
  gen_idx,
  rep_idx,
  col_idx = NULL,
  main_idx = NULL,
  design_type = c("RCBD", "LSD", "SPD"),
  method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"),
  tolerance = 1e-06
)
```

## Arguments

- data_mat:

  Numeric matrix with observations (rows) by traits (columns). May
  contain missing values (NA, NaN, Inf).

- gen_idx:

  Integer vector indicating genotype/treatment index for each
  observation (sub-plot in SPD).

- rep_idx:

  Integer vector indicating replicate/block (RCBD) or row (LSD) index
  for each observation.

- col_idx:

  Integer vector indicating column index for each observation (required
  for LSD only).

- main_idx:

  Integer vector indicating main plot treatment index for each
  observation (required for SPD only).

- design_type:

  Character string specifying experimental design: "RCBD" (default),
  "LSD" (Latin Square), or "SPD" (Split Plot).

- method:

  Character string specifying the estimation method. One of:

  REML

  :   Restricted Maximum Likelihood with variance components and BLUP
      (Best Linear Unbiased Prediction). Most robust for complex missing
      patterns. Default method. Iterations: 100. (RCBD, LSD only)

  Yates

  :   Traditional iterative formula: (r\*T + t\*B - G) / ((r-1)(t-1)).
      Simple and fast. Good for simple missing patterns. Iterations: 50.
      (RCBD, LSD only)

  Healy

  :   Healy & Westmacott method using weighted adjustment of treatment
      and block means. More stable than Yates for multiple missing
      values. Iterations: 50. (RCBD, LSD only)

  Regression

  :   Linear regression on treatment and block effects with QR
      decomposition. Non-iterative, single-pass estimation. Fast and
      stable. (RCBD, LSD only)

  Mean

  :   Simple mean substitution using treatment and block effects.
      Non-iterative, fastest method. Good for quick estimation. (RCBD,
      LSD, SPD)

  Bartlett

  :   ANCOVA using other traits as covariates with QR decomposition.
      Best when traits are correlated. Iterations: 30. (RCBD, LSD only)

- tolerance:

  Numeric convergence criterion. Iteration stops when maximum change in
  estimated values is below this threshold. Default: 1e-6.

## Value

A numeric matrix of the same dimensions as `data_mat` with all missing
values replaced by estimates.

## Details

The function handles missing values in RCBD, LSD, or SPD experiments by
iteratively estimating them until convergence. For RCBD, uses 2-way
blocking (genotypes × blocks); for LSD, uses 3-way blocking (genotypes ×
rows × columns); for SPD, uses nested structure (blocks \> main plots \>
sub-plots). Each method has different strengths:

**REML Method:** Uses variance component estimation with restricted
maximum likelihood. Applies shrinkage to predictions based on estimated
variance ratios. Most computationally intensive but handles complex
patterns well.

**Yates Method:** Uses the classical formula accounting for treatment,
block, and grand totals. Fast and simple, suitable for balanced designs
with few missing values.

**Healy & Westmacott Method:** Uses weighted adjustment based on
residuals from treatment and block effects. More stable than Yates when
multiple values are missing from same treatment or block. Particularly
effective for unbalanced missing patterns.

**Regression Method:** Fits linear model using complete observations
with treatment and block as factors. Non-iterative single-pass
estimation using QR decomposition for numerical stability. Fast and
deterministic - no convergence iterations needed.

**Mean Substitution Method:** Replaces missing values with treatment
mean + block mean - grand mean. Non-iterative, fastest method. Simple
additive model without regression. Best for quick estimation when
precision is less critical.

**Bartlett Method:** Performs ANCOVA using other traits as covariates to
predict missing values. Leverages trait correlations for better
predictions when traits are related. Uses QR decomposition for numerical
stability.

**Method Availability by Design:**

- RCBD: REML, Yates, Healy, Regression, Mean, Bartlett

- LSD: REML, Yates, Healy, Regression, Mean, Bartlett

- SPD: Mean only (other methods fall back to Mean)

## References

Yates, F. (1933). The analysis of replicated experiments when the field
results are incomplete. *Empire Journal of Experimental Agriculture*, 1,
129-142.

Healy, M. J. R., & Westmacott, M. (1956). Missing values in experiments
analysed on automatic computers. *Applied Statistics*, 5(3), 203-206.

Bartlett, M. S. (1937). Some examples of statistical methods of research
in agriculture and applied biology. *Supplement to the Journal of the
Royal Statistical Society*, 4(2), 137-183.
