# Comprehensive Example: Stochastic Simulation of Selection Indices
# This example demonstrates the use of Chapter 10 simulation features

# Load package
library(devtools)
load_all()

# ==============================================================================
# Example 1: Haldane's Mapping Function
# ==============================================================================
cat("Example 1: Haldane's Mapping Function\n")
cat("======================================\n\n")

# Convert genetic distances to recombination fractions
distances <- c(0, 0.05, 0.1, 0.2, 0.5, 1.0)
recomb_fractions <- haldane_mapping(distances)

cat("Genetic Distance (Morgans) -> Recombination Fraction\n")
for (i in seq_along(distances)) {
  cat(sprintf("  %.2f Morgan -> r = %.4f\n", distances[i], recomb_fractions[i]))
}

cat("\nNote: As distance increases, recombination approaches 0.5 (independent assortment)\n\n")


# ==============================================================================
# Example 2: Inverse Haldane Mapping
# ==============================================================================
cat("Example 2: Inverse Haldane Mapping\n")
cat("==================================\n\n")

# Convert recombination fractions back to distances
r_values <- c(0, 0.1, 0.2, 0.3, 0.4, 0.45)
d_values <- inverse_haldane_mapping(r_values)

cat("Recombination Fraction -> Genetic Distance (Morgans)\n")
for (i in seq_along(r_values)) {
  cat(sprintf("  r = %.2f -> %.4f Morgan\n", r_values[i], d_values[i]))
}
cat("\n\n")


# ==============================================================================
# Example 3: Small-Scale Simulation (10 cycles)
# ==============================================================================
cat("Example 3: Small-Scale Simulation\n")
cat("==================================\n\n")

set.seed(42)
results_small <- simulate_selection_cycles(
  n_cycles = 10,
  n_individuals = 200,
  n_loci = 30,
  n_traits = 3,
  heritability = 0.6,
  selection_proportion = 0.15,
  economic_weights = c(10, 5, 3),
  seed = 42
)

cat("Simulation completed!\n\n")
print(results_small)
cat("\n")
summary(results_small)
cat("\n\n")


# ==============================================================================
# Example 4: Simulation with Restricted Traits
# ==============================================================================
cat("Example 4: Restricted Selection (Trait 2 held at zero gain)\n")
cat("===========================================================\n\n")

set.seed(123)
results_restricted <- simulate_selection_cycles(
  n_cycles = 10,
  n_individuals = 200,
  n_loci = 30,
  n_traits = 3,
  heritability = 0.6,
  selection_proportion = 0.15,
  economic_weights = c(10, 5, 8),
  restricted_traits = 2,  # Restrict trait 2 to zero gain
  seed = 123
)

cat("Note: RLPSI and RESIM should show minimal gain for Trait 2\n\n")
summary(results_restricted)
cat("\n\n")


# ==============================================================================
# Example 5: Long-Term Simulation (Chapter 10: 50 cycles)
# ==============================================================================
cat("Example 5: Long-Term Simulation (50 cycles, as in Chapter 10)\n")
cat("=============================================================\n\n")
cat("This may take several minutes...\n\n")

set.seed(999)
results_longterm <- simulate_selection_cycles(
  n_cycles = 50,
  n_individuals = 1000,
  n_loci = 100,
  n_traits = 4,
  heritability = 0.5,
  selection_proportion = 0.1,
  economic_weights = c(12, 8, 6, 4),
  genetic_distances = rep(0.1, 99),  # 10 cM spacing
  seed = 999
)

cat("\nLong-term simulation completed!\n\n")
cat("Final Results (Cycle 50):\n")
final_cycle <- 50

cat("\nMean Genetic Values at Cycle 50:\n")
cat(sprintf("  LPSI:  %s\n", paste(round(results_longterm$lpsi_mean[final_cycle, ], 2), collapse = ", ")))
cat(sprintf("  ESIM:  %s\n", paste(round(results_longterm$esim_mean[final_cycle, ], 2), collapse = ", ")))
cat(sprintf("  RLPSI: %s\n", paste(round(results_longterm$rlpsi_mean[final_cycle, ], 2), collapse = ", ")))
cat(sprintf("  RESIM: %s\n", paste(round(results_longterm$resim_mean[final_cycle, ], 2), collapse = ", ")))

cat("\n\nTotal Cumulative Gain (50 cycles):\n")
cat(sprintf("  LPSI:  %s\n", paste(round(colSums(results_longterm$lpsi_gain), 2), collapse = ", ")))
cat(sprintf("  ESIM:  %s\n", paste(round(colSums(results_longterm$esim_gain), 2), collapse = ", ")))
cat(sprintf("  RLPSI: %s\n", paste(round(colSums(results_longterm$rlpsi_gain), 2), collapse = ", ")))
cat(sprintf("  RESIM: %s\n", paste(round(colSums(results_longterm$resim_gain), 2), collapse = ", ")))

cat("\n\n==============================================================\n")
cat("All examples completed successfully!\n")
cat("==============================================================\n")
