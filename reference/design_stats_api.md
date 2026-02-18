# Design Statistics API - Single Engine for ANOVA Computations

Unified API for experimental design statistics and ANOVA computations.
Replaces ad-hoc ANOVA calculations throughout the package, providing a
single source of truth for correction factors, sums of products, degrees
of freedom, and mean squares.

This function is a wrapper around design_stats() that computes
multivariate ANOVA statistics (mean square matrices) for all trait pairs
simultaneously.

## Usage

``` r
design_stats_api(
  data_mat,
  gen_idx,
  rep_idx,
  col_idx = NULL,
  main_idx = NULL,
  design_type = 1L
)
```

## Arguments

- data_mat:

  Numeric matrix of trait data (n_obs x n_traits)

- gen_idx:

  Integer vector of genotype indices (sub-plot treatments in SPD)

- rep_idx:

  Integer vector of replication/block indices (RCBD) or row indices
  (LSD)

- col_idx:

  Integer vector of column indices (for LSD, optional)

- main_idx:

  Integer vector of main plot indices (for SPD, optional)

- design_type:

  Integer design code: 1=RCBD, 2=LSD, 3=SPD

## Value

List with components compatible with legacy .calculate_anova():

- GMS:

  Genotype mean squares vector (diagonal of MSG)

- EMS:

  Error mean squares vector (diagonal of MSE)

- EMS_MAIN:

  Main plot error mean squares vector (SPD only, diagonal of MSG_MAIN)

- DFG:

  Degrees of freedom for genotypes/sub-plots

- DFE:

  Degrees of freedom for error (sub-plot error for SPD)

- DFE_MAIN:

  Degrees of freedom for main plot error (SPD only)

- n_rep:

  Number of replications

- n_gen:

  Number of genotypes/sub-plot treatments

- n_main:

  Number of main plot treatments (SPD only)

- MSG:

  Genotype mean square matrix (n_traits x n_traits)

- MSE:

  Error mean square matrix (n_traits x n_traits)

## Details

This function centralizes ANOVA computation logic, eliminating code
duplication across varcov.R, mean_performance.R, and other modules.

\*\*Design-specific formulas:\*\*

RCBD: - MSG = GSP / (g - 1) - MSE = ESP / ((g - 1) \* (r - 1))

LSD: - MSG = GSP / (t - 1) - MSE = ESP / ((t - 1) \* (t - 2))

SPD: - MSG = GSP / (b - 1) \[sub-plot treatments\] - MSE = ESP / (a \*
(b - 1) \* r) \[sub-plot error\] - MSE_MAIN = ESP_MAIN / ((a - 1) \*
(r - 1)) \[main plot error\]
