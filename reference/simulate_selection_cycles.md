# Simulate Multi-Cycle Selection Using Different Indices

Performs a stochastic simulation comparing the long-term performance of
four selection indices (LPSI, ESIM, RLPSI, RESIM) over multiple breeding
cycles. Models linkage between loci using Haldane's mapping function.

## Usage

``` r
simulate_selection_cycles(
  n_cycles = 50,
  n_individuals = 1000,
  n_loci = 100,
  n_traits = 3,
  heritability = 0.5,
  selection_proportion = 0.1,
  economic_weights = NULL,
  restricted_traits = NULL,
  genetic_distances = NULL,
  qtl_effects = NULL,
  seed = NULL
)
```

## Arguments

- n_cycles:

  Number of selection cycles to simulate (default: 50)

- n_individuals:

  Initial population size (default: 1000)

- n_loci:

  Number of QTL per trait (default: 100)

- n_traits:

  Number of quantitative traits (default: 3)

- heritability:

  Heritability for all traits (scalar or vector, default: 0.5)

- selection_proportion:

  Proportion of individuals selected (default: 0.1)

- economic_weights:

  Economic weights for LPSI (default: equal weights)

- restricted_traits:

  Trait indices to restrict to zero gain for RLPSI/RESIM (default: NULL)

- genetic_distances:

  Vector of genetic distances between loci in Morgans (default: 0.1)

- qtl_effects:

  Optional matrix of QTL effects (n_loci x n_traits)

- seed:

  Random seed for reproducibility (default: NULL)

## Value

List with components:

- `lpsi_gain` - Matrix of genetic gains per cycle for LPSI (n_cycles x
  n_traits)

- `esim_gain` - Matrix of genetic gains per cycle for ESIM

- `rlpsi_gain` - Matrix of genetic gains per cycle for RLPSI

- `resim_gain` - Matrix of genetic gains per cycle for RESIM

- `lpsi_mean` - Mean genetic value per cycle for LPSI (n_cycles x
  n_traits)

- `esim_mean` - Mean genetic value per cycle for ESIM

- `rlpsi_mean` - Mean genetic value per cycle for RLPSI

- `resim_mean` - Mean genetic value per cycle for RESIM

- `parameters` - List of simulation parameters used

## Details

**Simulation Procedure (Chapter 10):**

1\. Initialize population with random QTL alleles at linked loci 2. For
each cycle: - Compute genetic values from diploid genotypes - Add
environmental noise to create phenotypes - Calculate variance-covariance
matrices - Apply each selection index (LPSI, ESIM, RLPSI, RESIM) -
Select top individuals based on index scores - Generate offspring
through recombination (using Haldane's function) 3. Track genetic gain
and population mean across cycles

**The Four Indices (Chapter 10, Section 10.2):**

- **LPSI**: Unrestricted index maximizing genetic gain: b = P^(-1)Gw

- **ESIM**: Unrestricted eigen index maximizing accuracy: (P^(-1)C -
  lambda^2 I)b = 0

- **RLPSI**: Restricted LPSI with constraints: U'Gb = 0

- **RESIM**: Restricted ESIM with constraints: U'Cb = 0

**Important Simulation Assumptions:**

**Note:** Environmental variance is calculated once at generation 0 and
held constant across all cycles. Heritability will naturally decline
over cycles as genetic variance is depleted (Bulmer effect). This models
the biological reality that selection exhausts additive genetic variance
while environmental variation remains stable.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic simulation with 3 traits over 20 cycles
results <- simulate_selection_cycles(
  n_cycles = 20,
  n_individuals = 500,
  n_loci = 50,
  n_traits = 3,
  heritability = 0.5,
  selection_proportion = 0.1,
  economic_weights = c(10, 5, 3),
  seed = 123
)

# Plot genetic gains
plot(1:20, results$lpsi_mean[, 1],
  type = "l", col = "blue",
  ylab = "Mean Genetic Value", xlab = "Cycle",
  main = "Genetic Gain - Trait 1"
)
lines(1:20, results$esim_mean[, 1], col = "red")
lines(1:20, results$rlpsi_mean[, 1], col = "green")
lines(1:20, results$resim_mean[, 1], col = "orange")
legend("topleft", c("LPSI", "ESIM", "RLPSI", "RESIM"),
  col = c("blue", "red", "green", "orange"), lty = 1
)

# Restrict trait 2 to zero gain
results_restricted <- simulate_selection_cycles(
  n_cycles = 20,
  n_traits = 3,
  restricted_traits = 2,
  seed = 456
)
} # }
```
