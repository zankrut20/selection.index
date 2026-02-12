# Experimental Design Statistics Engine

Modular engine for experimental design statistics supporting RCBD, Latin
Square, and Split Plot designs. Computes correction factors, sums of
products, mean products, and degrees of freedom for variance-covariance
analysis and ANOVA statistics.

This function eliminates code duplication across gen.varcov(),
phen.varcov(), and mean.performance() by providing a centralized,
optimized implementation.

## Usage

``` r
design_stats(
  trait1,
  trait2 = trait1,
  genotypes,
  replications,
  columns = NULL,
  main_plots = NULL,
  design_type = c("RCBD", "LSD", "SPD"),
  calc_type = c("all", "sums_of_products", "mean_products", "anova_stats")
)
```

## Arguments

- trait1:

  Numeric vector of first trait observations

- trait2:

  Numeric vector of second trait observations (default: trait1 for
  variance)

- genotypes:

  Integer vector of genotype/treatment indices (sub-plot treatments in
  SPD)

- replications:

  Integer vector of replication/block indices (for RCBD) or row indices
  (for LSD)

- columns:

  Integer vector of column indices (required for Latin Square Design
  only)

- main_plots:

  Integer vector of main plot treatment indices (required for Split Plot
  Design only)

- design_type:

  Character string specifying design type: "RCBD" (default), "LSD"
  (Latin Square), or "SPD" (Split Plot)

- calc_type:

  Character string specifying calculation type:

  - `"sums_of_products"` - Returns CF, TSP, GSP, RSP, ESP (for
    covariance)

  - `"mean_products"` - Returns GMP, EMP (for variance components)

  - `"all"` - Returns all components (default)

  - `"anova_stats"` - Returns degrees of freedom and mean squares

## Value

List containing design statistics based on calc_type and design:

- `CF` - Correction factor

- `TSP` - Total sum of products

- `GSP` - Genotype/Sub-plot sum of products

- `RSP` - Replication sum of products

- `MSP` - Main plot sum of products (SPD only)

- `IMSP` - Main plot × Replication interaction SP (SPD only)

- `ESP` - Error sum of products (sub-plot error for SPD)

- `ESP_MAIN` - Main plot error sum of products (SPD only)

- `GMP` - Genotype mean product

- `EMP` - Error mean product (sub-plot error for SPD)

- `EMP_MAIN` - Main plot error mean product (SPD only)

- `DFG` - Degrees of freedom for genotypes/sub-plots

- `DFR` - Degrees of freedom for replications

- `DFM` - Degrees of freedom for main plots (SPD only)

- `DFE` - Degrees of freedom for error (sub-plot error for SPD)

- `DFE_MAIN` - Degrees of freedom for main plot error (SPD only)

- `n_genotypes` - Number of genotypes/sub-plot treatments

- `n_replications` - Number of replications

- `n_main_plots` - Number of main plot treatments (SPD only)

## Details

The function uses optimized Base R operations:

- [`rowsum()`](https://rdrr.io/r/base/rowsum.html) for fast grouped
  summations (5-10x faster than tapply)

- [`crossprod()`](https://rdrr.io/r/base/crossprod.html) for efficient
  matrix products (faster than sum(x\*y))

- Pre-computed constants to avoid repeated calculations

\*\*RCBD Model:\*\* Y_ij = μ + τ_i + β_j + ε_ij where τ_i = genotype
effect, β_j = block effect, ε_ij = error

\*\*LSD Model:\*\* Y_ijk = μ + τ_i + ρ_j + γ_k + ε_ijk where τ_i =
genotype effect, ρ_j = row effect, γ_k = column effect, ε_ijk = error

\*\*SPD Model:\*\* Y_ijk = μ + ρ_i + α_j + δ_ij + τ_k + (ατ)\_jk + ε_ijk
where ρ_i = block effect, α_j = main plot effect, δ_ij = main plot
error, τ_k = sub-plot effect (genotype), (ατ)\_jk = interaction, ε_ijk =
sub-plot error

## References

Cochran, W. G., & Cox, G. M. (1957). Experimental designs (2nd ed.).
Wiley.

Steel, R. G. D., & Torrie, J. H. (1980). Principles and procedures of
statistics: A biometrical approach (2nd ed.). McGraw-Hill.

Gomez, K. A., & Gomez, A. A. (1984). Statistical procedures for
agricultural research (2nd ed.). Wiley.
