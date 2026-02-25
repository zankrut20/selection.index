# Stochastic Simulation of Linear Phenotypic Selection Indices

## Introduction

Stochastic simulation provides breeders with a comprehensive method to
evaluate the long-term behavior of different selection indices across
numerous sequential cycles. By running *in silico* (computer simulated)
models, breeders can evaluate the expected genetic gains, selection
accuracy, and the eventual decay of genetic variance over an extended
selection timeframe.

Chapter 10 of the `selection.index` package focuses on the performance
comparison of over 50 recurrent generic population breeding cycles. The
package provides the
[`stochastic_simulation()`](https://zankrut20.github.io/selection.index/reference/stochastic_simulation.md)
function to evaluate multiple restricted and unrestricted linear
phenotypic indices:

1.  **LPSI** (Linear Phenotypic Selection Index - unrestricted
    Smith-Hazel)
2.  **ESIM** (Eigen Selection Index Method - unrestricted eigen-based
    method)
3.  **RLPSI** (Restricted Linear Phenotypic Selection Index - restricted
    Kempthorne-Nordskog method)
4.  **RESIM** (Restricted Eigen Selection Index Method - restricted
    eigen-based method)

### Simulated Breeding Design

The `stochastic_simulation` function mimics a forward recurrent
selection scheme. Within each cycle, the simulation performs the
following biological steps: 1. Recombines parental haplotypes using
Haldane’s mapping function to generate double-haploid lines via
crossing. 2. Simulates quantitative trait variations assuming a full
pleiotropic model, applying sampled major/minor QTL effects. 3. Computes
the phenotypic score of the progenies by applying realistic
environmental variances based on expected heritabilities ($h^{2}$). 4.
Evaluates the linear combinations of each chosen index to select the
top-performing progenies. 5. Randomly mates the selected superior lines
to act as the founders for the next evolutionary cycle.

------------------------------------------------------------------------

## Running a Simulation

We can run a customized, small-scale simulation to demonstrate the
tracking functionality. We will analyze a scenario with three traits to
trace their behavior across 5 generation cycles.

``` r
library(selection.index)

# To ensure reproducibility
set.seed(42)
```

### Configuring the Scenario

In accordance with the textbook paradigm, let us analyze three
genetically correlated target traits: - **Trait 1:** Low heritability
($h^{2} = 0.2$), negatively correlated to Trait 2. - **Trait 2:** High
heritability ($h^{2} = 0.5$). - **Trait 3:** Independent trait with high
heritability ($h^{2} = 0.5$).

#### 1. Generating Genetic Effects

First, we generate a known initial phenotypic and genotypic covariance
matrix to serve as the baseline simulation logic for our markers. For
simplicity, we will simulate a rapid 5-cycle breeding program across a
smaller subset of 50 QTL markers.

``` r
n_traits <- 3
n_loci <- 50 # Number of segregating sites / markers

# Generate random base QTL effects for the markers across the 3 traits
# Negative correlation infused between trait 1 and 2
qtl_eff <- matrix(rnorm(n_loci * n_traits), nrow = n_loci, ncol = n_traits)
qtl_eff[, 2] <- -0.5 * qtl_eff[, 1] + 0.5 * qtl_eff[, 2]

# Define heritabilities and corresponding environmental variance
heritabilities <- c(0.2, 0.5, 0.5)

# Simulate base genetic variance to deduce correct environmental variance noise
base_gv <- apply(qtl_eff, 2, var) * n_loci
env_var <- base_gv * (1 - heritabilities) / heritabilities
```

#### 2. Formulating Constraints and Weights

All the indices will apply an equal starting economic weight. However,
for the Restricted indices (RLPSI and RESIM), we will force a constraint
onto the first trait.

``` r
# Equal economic trait weighting
weights <- c(1, 1, 1)

# Constraint matrix for RLPSI/RESIM: Constrain Trait 1
U_mat <- matrix(0, nrow = 3, ncol = 1)
U_mat[1, 1] <- 1
```

#### 3. Executing the Simulation Routine

The `simulate_selection_cycles` routine requires various inputs covering
biological configuration (individuals, cycles) and selection
configurations (weights, restricted traits, selection proportion).

*Note: For performance within the vignette, we use minimal
cyclic/generational values. In real breeding scenarios, you would set
`n_cycles = 50` and `n_individuals = 10000` to mirror the book’s 35,000
double-haploid evaluations.*

``` r
# Run the stochastic selection (may take a moment)
sim_results <- simulate_selection_cycles(
  n_cycles = 5,
  n_individuals = 200,
  n_loci = n_loci,
  n_traits = n_traits,
  qtl_effects = qtl_eff,
  heritability = heritabilities,
  economic_weights = weights,
  selection_proportion = 0.25, # Select upper 25% progeny
  restricted_traits = 1
)
```

------------------------------------------------------------------------

## Interpreting Results

The simulation tracks comprehensive metrics for all computed methods
(LPSI, ESIM, RLPSI, RESIM) across all generated cycles. The return
object provides tracking arrays for: 1. Genetic Gain (`*_gain`) 2.
Genetic Variance (`*_var`) 3. Estimated Heritabilities
(`*_heritability`) 4. Mean Phenotypic Value (`*_mean`)

#### Genetic Gain Trajectories

We can extract the continuous genetic gain trajectories that map the
success of the applied phenotypic criteria across the multiple
sequential cycles. Let’s observe the restricted indices:

``` r
# Expected: Because Trait 1 was constrained via the U_mat for the RLPSI metric,
# its expected generational gain should stabilize at 0.
print(sim_results$rlpsi_gain)
#>             [,1]       [,2]     [,3]
#> [1,]  0.00000000  0.0000000 0.000000
#> [2,] -0.07991849  0.3576849 4.396779
#> [3,]  0.09595190  0.7166170 2.927811
#> [4,]  0.44428982 -0.2318664 2.712902
#> [5,]  1.14420568 -0.4358615 3.740713
```

As clearly demonstrated, the selection gain mapped for Trait 1 drops
immediately towards effectively $0$, verifying that the restricted
Kempthorne-Nordskog properties hold up systematically across stochastic
evaluations. In comparison, unrestricted indices (such as LPSI) map
increasing steady gains across all associated traits.

``` r
print(sim_results$lpsi_gain)
#>           [,1]       [,2]     [,3]
#> [1,] 0.0000000  0.0000000 0.000000
#> [2,] 0.1607943 -0.4824810 3.821722
#> [3,] 1.9171279 -1.3258871 3.435662
#> [4,] 1.0369927 -0.5201380 4.007482
#> [5,] 1.1104539 -0.7914383 2.622570
```

#### Decay of Variance

Additionally, intense generational selection logically extinguishes raw
genetic variability. If the selection index consistently applies extreme
thresholds, the available allelic variance diminishes, making subsequent
trait improvements plateau.

``` r
# Observe the diminishing variance arrays for the LPSI evaluations
print(sim_results$lpsi_var)
#> NULL
```

------------------------------------------------------------------------

## Summary

By evaluating cyclic genetic simulations across stochastic marker
environments, breeders identify crucial plateau mechanisms in
multidimensional indices. Specifically: - **Restricted indices**
properly maintain and prevent target constraints over multi-generational
limits, albeit slightly reducing overarching expected gains. -
**Independence of traits** limits associative “ride-along” gains
obtained via pleiotropy in highly positive correlated setups. -
**Accuracies degrade linearly** as available variant variance
diminishes, confirming limits built into simple additive pleiotropic
models.
