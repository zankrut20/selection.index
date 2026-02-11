# =============================================================================
# PHASE 1 DEMONSTRATION: Enhanced Desired Gains Index
# =============================================================================
# This script demonstrates the new features added to dg_lpsi() in Phase 1:
# 1. Implied economic weights calculation (ŵ = G^-1 P b)
# 2. Feasibility checking for desired gains
# 3. Enhanced print and summary methods
# =============================================================================

library(selection.index)

# Load example data
data(seldata)

# Calculate variance-covariance matrices
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# =============================================================================
# EXAMPLE 1: Basic usage with realistic gains
# =============================================================================
cat("\n=== EXAMPLE 1: Realistic Desired Gains ===\n")

# Set realistic desired gains (small values)
desired_gains <- c(0.5, 0.3, 0.4, 0.2, 0.3, 0.5, 0.4)

result1 <- dg_lpsi(
  pmat = pmat,
  gmat = gmat,
  d = desired_gains,
  return_implied_weights = TRUE,
  check_feasibility = TRUE
)

# Print the result (compact view)
print(result1)

# Get detailed summary
cat("\n--- Detailed Summary ---\n")
summary(result1)

# =============================================================================
# EXAMPLE 2: Unrealistic gains trigger warnings
# =============================================================================
cat("\n\n=== EXAMPLE 2: Unrealistic Desired Gains ===\n")

# Set unrealistically large desired gains
unrealistic_gains <- c(5, 6, 7, 4, 5, 8, 6)

result2 <- dg_lpsi(
  pmat = pmat,
  gmat = gmat,
  d = unrealistic_gains,
  return_implied_weights = TRUE,
  check_feasibility = TRUE,
  selection_intensity = 2.063  # Standard selection intensity (10% selected)
)

# Note: This will produce warnings about unrealistic gains

# Examine feasibility metrics
cat("\n--- Feasibility Metrics ---\n")
print(result2$feasibility)

# =============================================================================
# EXAMPLE 3: Implied economic weights interpretation
# =============================================================================
cat("\n\n=== EXAMPLE 3: Understanding Implied Weights ===\n")

# Set selective gains (favor some traits over others)
selective_gains <- c(1.0, 0.1, 0.2, 0.1, 0.5, 0.8, 0.3)

result3 <- dg_lpsi(
  pmat = pmat,
  gmat = gmat,
  d = selective_gains,
  return_implied_weights = TRUE,
  check_feasibility = FALSE  # Disable feasibility checking
)

cat("\n--- Implied Economic Weights ---\n")
cat("These weights show the 'cost' of improving each trait.\n")
cat("Large absolute values indicate traits that are expensive to improve.\n\n")

weights_comparison <- data.frame(
  Trait = colnames(gmat),
  Desired_Gain = selective_gains,
  Implied_Weight = round(result3$implied_weights, 3),
  Normalized_Weight = round(result3$implied_weights_normalized, 3)
)
print(weights_comparison)

# =============================================================================
# EXAMPLE 4: High selection intensity scenario
# =============================================================================
cat("\n\n=== EXAMPLE 4: High Selection Intensity ===\n")

# Use higher selection intensity (top 1% selected, i ≈ 2.67)
moderate_gains <- c(2, 2.5, 3, 1.5, 2, 3, 2.5)

result4 <- dg_lpsi(
  pmat = pmat,
  gmat = gmat,
  d = moderate_gains,
  selection_intensity = 2.67,  # Top 1% of population
  check_feasibility = TRUE,
  return_implied_weights = TRUE
)

cat("\n--- Feasibility with Intense Selection ---\n")
cat("With higher selection intensity, larger gains become feasible.\n\n")
print(result4$feasibility[, c("trait", "feasibility_ratio", "is_realistic")])

# =============================================================================
# EXAMPLE 5: Comparing different gain scenarios
# =============================================================================
cat("\n\n=== EXAMPLE 5: Comparing Scenarios ===\n")

# Scenario A: Equal gains for all traits
equal_gains <- rep(0.5, ncol(gmat))
result_equal <- dg_lpsi(pmat, gmat, equal_gains, check_feasibility = FALSE)

# Scenario B: Prioritize yield traits
yield_priority <- c(1.5, 0.2, 0.3, 0.2, 0.5, 1.0, 0.3)
result_yield <- dg_lpsi(pmat, gmat, yield_priority, check_feasibility = FALSE)

cat("\n--- Index Efficiency Comparison ---\n")
comparison <- data.frame(
  Scenario = c("Equal Gains", "Yield Priority"),
  Heritability = c(result_equal$hI2, result_yield$hI2),
  Correlation_HI = c(result_equal$rHI, result_yield$rHI),
  Expected_Genetic_Gain = c(result_equal$Delta_G, result_yield$Delta_G)
)
print(comparison)

cat("\n--- Implied Weight Tradeoffs ---\n")
cat("\nEqual Gains:\n")
print(round(result_equal$implied_weights_normalized, 3))
cat("\nYield Priority:\n")
print(round(result_yield$implied_weights_normalized, 3))

# =============================================================================
# KEY INSIGHTS FROM PHASE 1 ENHANCEMENTS
# =============================================================================
cat("\n\n=== KEY INSIGHTS ===\n")
cat("1. Implied weights reveal the genetic/phenotypic cost structure\n")
cat("2. Feasibility checking prevents unrealistic breeding targets\n")
cat("3. Selection intensity affects what gains are achievable\n")
cat("4. Different gain priorities result in different weight tradeoffs\n")
cat("5. Enhanced output helps interpret and validate breeding strategies\n\n")
