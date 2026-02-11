# =============================================================================
# PHASE 2 DEMONSTRATION: Base Index (Williams, 1962)
# =============================================================================
# This script demonstrates the Base Index implementation where coefficients
# equal economic weights (b = w). This is a simple, unoptimized approach
# useful when covariance matrices are unreliable or as a baseline comparison.
# =============================================================================

library(selection.index)

# Load example data
data(seldata)
data(weight)

# Calculate variance-covariance matrices
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# =============================================================================
# EXAMPLE 1: Basic Base Index
# =============================================================================
cat("\n=== EXAMPLE 1: Basic Base Index Usage ===\n")

# Define economic weights (equal weights for simplicity)
weights <- rep(1, ncol(gmat))

result1 <- base_index(
  pmat = pmat,
  gmat = gmat,
  wmat = weights,
  selection_intensity = 2.063,
  compare_to_lpsi = TRUE
)

# Print the result
print(result1)

# Get detailed summary
cat("\n--- Detailed Summary ---\n")
summary(result1)

# =============================================================================
# EXAMPLE 2: Comparison with Optimal LPSI
# =============================================================================
cat("\n\n=== EXAMPLE 2: Base Index vs Optimal LPSI ===\n")

# Define economic weights favoring some traits
weights_selective <- c(10, 5, 3, 3, 5, 8, 4)

result2 <- base_index(
  pmat = pmat,
  gmat = gmat,
  wmat = weights_selective,
  compare_to_lpsi = TRUE
)

# Compare efficiencies
cat("\n--- Efficiency Comparison ---\n")
cat(sprintf("Base Index GA:       %.4f\n", result2$GA))
cat(sprintf("Optimal LPSI GA:     %.4f\n", result2$lpsi_comparison$GA_lpsi))
cat(sprintf("Efficiency Ratio:    %.2f%% of LPSI\n", 
            result2$lpsi_comparison$efficiency_ratio * 100))

# Response comparison
cat("\n--- Response Comparison ---\n")
response_comparison <- data.frame(
  Trait = colnames(gmat),
  Base_DeltaG = round(result2$Delta_G, 4),
  LPSI_DeltaG = round(result2$lpsi_comparison$Delta_G_lpsi, 4),
  Difference = round(result2$lpsi_comparison$Delta_G_lpsi - result2$Delta_G, 4)
)
print(response_comparison)

# =============================================================================
# EXAMPLE 3: Multiple Weight Scenarios
# =============================================================================
cat("\n\n=== EXAMPLE 3: Multiple Weight Scenarios ===\n")

# Create weight matrix with different scenarios
weights_matrix <- cbind(
  Equal = rep(1, ncol(gmat)),              # Equal weights
  Yield_Focus = c(10, 2, 2, 1, 3, 8, 2),  # Focus on yield traits
  Quality_Focus = c(2, 8, 10, 8, 2, 3, 5) # Focus on quality traits
)

# Calculate Base Index for each scenario
results_equal <- base_index(pmat, gmat, weights_matrix, wcol = 1, compare_to_lpsi = FALSE)
results_yield <- base_index(pmat, gmat, weights_matrix, wcol = 2, compare_to_lpsi = FALSE)
results_quality <- base_index(pmat, gmat, weights_matrix, wcol = 3, compare_to_lpsi = FALSE)

# Compare genetic advances
cat("\n--- Scenario Comparison ---\n")
scenario_comparison <- data.frame(
  Scenario = c("Equal Weights", "Yield Focus", "Quality Focus"),
  GA = c(results_equal$GA, results_yield$GA, results_quality$GA),
  hI2 = c(results_equal$hI2, results_yield$hI2, results_quality$hI2),
  rHI = c(results_equal$rHI, results_yield$rHI, results_quality$rHI)
)
print(round(scenario_comparison, 4))

# =============================================================================
# EXAMPLE 4: Robustness to Unreliable Estimates
# =============================================================================
cat("\n\n=== EXAMPLE 4: Robustness Demonstration ===\n")

# Simulate scenario: perturb covariance matrices (simulate estimation error)
set.seed(123)
perturbation <- 0.1  # 10% perturbation

pmat_perturbed <- pmat * (1 + rnorm(length(pmat), 0, perturbation))
gmat_perturbed <- gmat * (1 + rnorm(length(gmat), 0, perturbation))

# Make sure they're still symmetric and positive definite
pmat_perturbed <- (pmat_perturbed + t(pmat_perturbed)) / 2
gmat_perturbed <- (gmat_perturbed + t(gmat_perturbed)) / 2

# Ensure positive definiteness (add small diagonal if needed)
diag(pmat_perturbed) <- diag(pmat_perturbed) + abs(min(eigen(pmat_perturbed)$values)) + 0.1
diag(gmat_perturbed) <- diag(gmat_perturbed) + abs(min(eigen(gmat_perturbed)$values)) + 0.1

weights <- rep(1, ncol(gmat))

# Original matrices
result_original <- base_index(pmat, gmat, weights, compare_to_lpsi = FALSE)

# Perturbed matrices
result_perturbed <- base_index(pmat_perturbed, gmat_perturbed, weights, compare_to_lpsi = FALSE)

cat("\n--- Stability Under Perturbation ---\n")
cat(sprintf("Original GA:     %.4f\n", result_original$GA))
cat(sprintf("Perturbed GA:    %.4f\n", result_perturbed$GA))
cat(sprintf("Relative change: %.2f%%\n", 
            abs(result_perturbed$GA - result_original$GA) / result_original$GA * 100))

cat("\nBase Index is robust because it doesn't require matrix inversion,\n")
cat("making it stable when covariance estimates are unreliable.\n")

# =============================================================================
# EXAMPLE 5: Selection Intensity Effects
# =============================================================================
cat("\n\n=== EXAMPLE 5: Selection Intensity Effects ===\n")

weights <- c(5, 3, 4, 2, 3, 6, 3)

# Different selection intensities
intensities <- c(1.40, 1.76, 2.06, 2.42, 2.67)  # 50%, 25%, 10%, 5%, 1% selected
names(intensities) <- c("50%", "25%", "10%", "5%", "1%")

intensity_results <- lapply(intensities, function(i) {
  base_index(pmat, gmat, weights, selection_intensity = i, compare_to_lpsi = FALSE)
})

cat("\n--- Selection Intensity Impact ---\n")
intensity_comparison <- data.frame(
  Selection_Prop = names(intensities),
  Intensity = intensities,
  GA = sapply(intensity_results, function(x) x$GA),
  Mean_DeltaG = sapply(intensity_results, function(x) mean(x$Delta_G))
)
print(round(intensity_comparison, 4))

cat("\nGenetic advance scales linearly with selection intensity.\n")
cat("More intense selection → larger genetic gains.\n")

# =============================================================================
# EXAMPLE 6: Using Weight Matrix from Package Data
# =============================================================================
cat("\n\n=== EXAMPLE 6: Using Package Weight Matrix ===\n")

# Use weight matrix from package data
result6 <- base_index(
  pmat = pmat,
  gmat = gmat,
  wmat = weight[,-1],  # Remove first column (trait names)
  wcol = 1,            # Use first weight column
  compare_to_lpsi = TRUE
)

print(result6)

# =============================================================================
# EXAMPLE 7: When Base Index Performs Well
# =============================================================================
cat("\n\n=== EXAMPLE 7: Scenarios Where Base Index Excels ===\n")

# Create scenario with weak correlations (diagonal-dominant)
P_uncorr <- diag(diag(pmat)) + 0.1 * (pmat - diag(diag(pmat)))
G_uncorr <- diag(diag(gmat)) + 0.1 * (gmat - diag(diag(gmat)))

weights <- rep(1, ncol(gmat))

result_uncorr <- base_index(P_uncorr, G_uncorr, weights, compare_to_lpsi = TRUE)

cat("\n--- Performance with Weak Correlations ---\n")
cat(sprintf("Efficiency vs LPSI: %.2f%%\n", 
            result_uncorr$lpsi_comparison$efficiency_ratio * 100))

if (result_uncorr$lpsi_comparison$efficiency_ratio > 0.90) {
  cat("\n✓ Base Index performs well (>90% of LPSI efficiency)\n")
  cat("  When correlations are weak, Base Index is nearly optimal.\n")
}

# =============================================================================
# KEY INSIGHTS FROM PHASE 2
# =============================================================================
cat("\n\n=== KEY INSIGHTS ===\n")
cat("1. Base Index (b = w) is simple: no matrix inversion required\n")
cat("2. Less efficient than LPSI but more robust to estimation errors\n")
cat("3. Performs well when correlations between traits are weak\n")
cat("4. Useful baseline for comparing optimized selection indices\n")
cat("5. Genetic advance scales linearly with selection intensity\n")
cat("6. Particularly valuable when sample sizes are small\n")
cat("7. Trade-off: simplicity and robustness vs optimal efficiency\n\n")
