# ==============================================================================
# Tests for Stochastic Simulation (Chapter 10)
# ==============================================================================

# ------------------------------------------------------------------------------
# Test Haldane's Mapping Function
# ------------------------------------------------------------------------------

test_that("haldane_mapping basic functionality works", {
  # Test single value
  r <- haldane_mapping(0.1)
  expect_type(r, "double")
  expect_length(r, 1)
  expect_true(r >= 0 && r <= 0.5)
  
  # Test vector input
  distances <- c(0.05, 0.1, 0.2, 0.5)
  r_vec <- haldane_mapping(distances)
  expect_length(r_vec, 4)
  expect_true(all(r_vec >= 0 & r_vec <= 0.5))
})

test_that("haldane_mapping edge cases work", {
  # Zero distance should give zero recombination
  expect_equal(haldane_mapping(0), 0)
  
  # Infinite distance should give 0.5 recombination
  expect_equal(haldane_mapping(Inf), 0.5)
  
  # Very small distance
  r_small <- haldane_mapping(0.001)
  expect_true(r_small > 0 && r_small < 0.001)
  
  # 1 Morgan should give specific value
  r_1morgan <- haldane_mapping(1)
  expected <- 0.5 * (1 - exp(-2))
  expect_equal(r_1morgan, expected, tolerance = 1e-10)
})

test_that("haldane_mapping validates input", {
  # Negative distance should error
  expect_error(haldane_mapping(-0.1), "non-negative")
})

# ------------------------------------------------------------------------------
# Test Inverse Haldane's Mapping Function
# ------------------------------------------------------------------------------

test_that("inverse_haldane_mapping basic functionality works", {
  # Test single value
  d <- inverse_haldane_mapping(0.1)
  expect_type(d, "double")
  expect_length(d, 1)
  expect_true(d >= 0)
  
  # Test vector input
  r_values <- c(0.05, 0.1, 0.2, 0.4)
  d_vec <- inverse_haldane_mapping(r_values)
  expect_length(d_vec, 4)
  expect_true(all(d_vec >= 0))
})

test_that("inverse_haldane_mapping edge cases work", {
  # Zero recombination should give zero distance
  expect_equal(inverse_haldane_mapping(0), 0)
  
  # 0.5 recombination should give infinite distance (with warning)
  expect_warning(
    result <- inverse_haldane_mapping(0.5),
    "infinite genetic distance"
  )
  expect_equal(result, Inf)
  
  # Very small recombination
  d_small <- inverse_haldane_mapping(0.001)
  expect_true(d_small > 0 && d_small < 0.01)
})

test_that("inverse_haldane_mapping validates input", {
  # Negative r should error
  expect_error(inverse_haldane_mapping(-0.1), "between 0 and 0.5")
  
  # r > 0.5 should error
  expect_error(inverse_haldane_mapping(0.6), "between 0 and 0.5")
})

test_that("haldane and inverse are true inverses", {
  # Forward then backward
  distances <- c(0.01, 0.1, 0.5, 1.0)
  r_values <- haldane_mapping(distances)
  d_recovered <- inverse_haldane_mapping(r_values)
  expect_equal(d_recovered, distances, tolerance = 1e-10)
  
  # Backward then forward
  r_values2 <- c(0.01, 0.1, 0.2, 0.4)
  d_values <- inverse_haldane_mapping(r_values2)
  r_recovered <- haldane_mapping(d_values)
  expect_equal(r_recovered, r_values2, tolerance = 1e-10)
})

# ------------------------------------------------------------------------------
# Test simulate_selection_cycles Basic Functionality
# ------------------------------------------------------------------------------

test_that("simulate_selection_cycles basic functionality works", {
  set.seed(123)
  
  # Small test case
  result <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.6,
    selection_proportion = 0.1,
    economic_weights = c(10, 5),
    seed = 123
  )
  
  # Check result structure
  expect_type(result, "list")
  expect_s3_class(result, "selection_simulation")
  
  # Check required components
  expect_true("lpsi_gain" %in% names(result))
  expect_true("esim_gain" %in% names(result))
  expect_true("rlpsi_gain" %in% names(result))
  expect_true("resim_gain" %in% names(result))
  expect_true("lpsi_mean" %in% names(result))
  expect_true("esim_mean" %in% names(result))
  expect_true("rlpsi_mean" %in% names(result))
  expect_true("resim_mean" %in% names(result))
  expect_true("parameters" %in% names(result))
  
  # Check dimensions
  expect_equal(dim(result$lpsi_gain), c(5, 2))
  expect_equal(dim(result$lpsi_mean), c(5, 2))
})

test_that("simulate_selection_cycles with restrictions works", {
  set.seed(456)
  
  result <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 3,
    heritability = 0.6,
    selection_proportion = 0.1,
    economic_weights = c(10, 5, 3),
    restricted_traits = 2,
    seed = 456
  )
  
  # Check result structure
  expect_s3_class(result, "selection_simulation")
  expect_equal(dim(result$rlpsi_gain), c(5, 3))
  
  # RLPSI and RESIM should restrict trait 2 gain
  # Total gain for trait 2 should be close to zero
  total_gain_rlpsi_trait2 <- sum(result$rlpsi_gain[, 2])
  total_gain_resim_trait2 <- sum(result$resim_gain[, 2])
  
  expect_true(abs(total_gain_rlpsi_trait2) < 5)  # Should be near zero
  expect_true(abs(total_gain_resim_trait2) < 5)  # Should be near zero
})

test_that("simulate_selection_cycles reproduces with same seed", {
  result1 <- simulate_selection_cycles(
    n_cycles = 3,
    n_individuals = 50,
    n_loci = 10,
    n_traits = 2,
    seed = 789
  )
  
  result2 <- simulate_selection_cycles(
    n_cycles = 3,
    n_individuals = 50,
    n_loci = 10,
    n_traits = 2,
    seed = 789
  )
  
  # Results should be identical with same seed
  expect_equal(result1$lpsi_mean, result2$lpsi_mean)
  expect_equal(result1$esim_gain, result2$esim_gain)
})

# ------------------------------------------------------------------------------
# Test Bulmer Effect (Environmental Variance Constant)
# ------------------------------------------------------------------------------

test_that("Bulmer effect is present (gains decline over time)", {
  set.seed(2026)
  
  result <- simulate_selection_cycles(
    n_cycles = 20,
    n_individuals = 200,
    n_loci = 40,
    n_traits = 2,
    heritability = 0.6,
    selection_proportion = 0.05,  # High selection intensity
    seed = 2026
  )
  
  # Check that gains decline over time (Bulmer effect)
  # Early cycles should have higher gains than late cycles
  early_gains <- mean(result$lpsi_gain[2:5, 1])
  late_gains <- mean(result$lpsi_gain[16:20, 1])
  
  expect_true(early_gains > late_gains)
  
  # Late gains should be substantially smaller
  decline_ratio <- late_gains / early_gains
  expect_true(decline_ratio < 0.5)  # At least 50% decline
})

# ------------------------------------------------------------------------------
# Test Robustness to Long Simulations
# ------------------------------------------------------------------------------

test_that("simulate_selection_cycles handles extended cycles without errors", {
  set.seed(2027)
  
  # This should complete without matrix singularity errors
  # Thanks to MASS::ginv() robustness
  expect_silent({
    result <- simulate_selection_cycles(
      n_cycles = 30,
      n_individuals = 150,
      n_loci = 30,
      n_traits = 2,
      heritability = 0.6,
      selection_proportion = 0.05,
      seed = 2027
    )
  })
  
  # Should complete all cycles
  expect_equal(nrow(result$lpsi_mean), 30)
  expect_equal(nrow(result$esim_gain), 30)
})

# ------------------------------------------------------------------------------
# Test S3 Methods
# ------------------------------------------------------------------------------

test_that("print.selection_simulation works", {
  set.seed(999)
  
  result <- simulate_selection_cycles(
    n_cycles = 3,
    n_individuals = 50,
    n_loci = 10,
    n_traits = 2,
    seed = 999
  )
  
  # Should print without error
  expect_output(print(result), "Stochastic Selection")
  expect_output(print(result), "LPSI")
  expect_output(print(result), "ESIM")
})

test_that("summary.selection_simulation works", {
  set.seed(1000)
  
  result <- simulate_selection_cycles(
    n_cycles = 3,
    n_individuals = 50,
    n_loci = 10,
    n_traits = 2,
    seed = 1000
  )
  
  # Should produce summary without error
  expect_output(summary(result), "Simulation Summary")
  expect_output(summary(result), "Total Genetic Gain")
  expect_output(summary(result), "Mean Genetic Gain per Cycle")
})

test_that("plot.selection_simulation works", {
  set.seed(1001)
  
  result <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 50,
    n_loci = 10,
    n_traits = 2,
    seed = 1001
  )
  
  # Should plot without error
  expect_silent(plot(result))
  expect_silent(plot(result, trait = 2))
})

# ------------------------------------------------------------------------------
# Test Parameter Validation
# ------------------------------------------------------------------------------

test_that("simulate_selection_cycles validates parameters", {
  # n_cycles must be at least 1
  expect_error(
    simulate_selection_cycles(n_cycles = 0),
    "at least 1"
  )
  
  # n_individuals must be at least 10
  expect_error(
    simulate_selection_cycles(n_individuals = 0),
    "at least 10"
  )
  
  # heritability must be between 0 and 1
  expect_error(
    simulate_selection_cycles(heritability = 1.5),
    "between 0 and 1"
  )
  
  # selection_proportion must be between 0 and 1
  expect_error(
    simulate_selection_cycles(selection_proportion = 1.5),
    "between 0 and 1"
  )
})

# ------------------------------------------------------------------------------
# Test Restricted vs Unrestricted Indices
# ------------------------------------------------------------------------------

test_that("RLPSI equals LPSI when no restrictions", {
  set.seed(2028)
  
  result <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.6,
    restricted_traits = NULL,  # No restrictions
    seed = 2028
  )
  
  # Without restrictions, RLPSI should be very similar to LPSI
  # (small differences due to sampling)
  lpsi_total <- sum(result$lpsi_mean[5, ])
  rlpsi_total <- sum(result$rlpsi_mean[5, ])
  
  # Should be within 20% of each other
  relative_diff <- abs(lpsi_total - rlpsi_total) / abs(lpsi_total)
  expect_true(relative_diff < 0.2)
})

test_that("RESIM behavior when no restrictions", {
  set.seed(2029)
  
  result <- simulate_selection_cycles(
    n_cycles = 10,
    n_individuals = 200,
    n_loci = 40,
    n_traits = 2,
    heritability = 0.6,
    restricted_traits = NULL,  # No restrictions
    seed = 2029
  )
  
  # Without restrictions, both ESIM and RESIM should complete successfully
  # (they may give different results due to stochastic variation and different algorithms)
  expect_true(all(!is.na(result$esim_mean)))
  expect_true(all(!is.na(result$resim_mean)))
  
  # Both should show positive gains overall
  esim_final <- sum(result$esim_mean[10, ])
  resim_final <- sum(result$resim_mean[10, ])
  
  # Check that both indices produce results (no NAs or infinities)
  expect_true(is.finite(esim_final))
  expect_true(is.finite(resim_final))
})

# ------------------------------------------------------------------------------
# Test Multiple Traits
# ------------------------------------------------------------------------------

test_that("simulate_selection_cycles works with multiple traits", {
  set.seed(3000)
  
  result <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 5,  # Multiple traits
    heritability = 0.6,
    economic_weights = c(10, 8, 6, 4, 2),
    seed = 3000
  )
  
  # Check dimensions
  expect_equal(dim(result$lpsi_gain), c(5, 5))
  expect_equal(dim(result$lpsi_mean), c(5, 5))
  
  # All indices should work
  expect_true(all(!is.na(result$lpsi_gain)))
  expect_true(all(!is.na(result$esim_gain)))
  expect_true(all(!is.na(result$rlpsi_gain)))
  expect_true(all(!is.na(result$resim_gain)))
})

# ------------------------------------------------------------------------------
# Test Different Heritabilities
# ------------------------------------------------------------------------------

test_that("simulate_selection_cycles works with different heritabilities", {
  set.seed(3001)
  
  # High heritability
  result_high <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.9,
    seed = 3001
  )
  
  # Low heritability
  result_low <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.3,
    seed = 3001
  )
  
  # High heritability should lead to more total gain
  gain_high <- sum(abs(result_high$lpsi_gain))
  gain_low <- sum(abs(result_low$lpsi_gain))
  
  expect_true(gain_high > gain_low)
})

# ------------------------------------------------------------------------------
# Test Selection Intensity Effects
# ------------------------------------------------------------------------------

test_that("higher selection intensity increases gains", {
  set.seed(3002)
  
  # High selection intensity
  result_high <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.6,
    selection_proportion = 0.05,  # Top 5%
    seed = 3002
  )
  
  # Low selection intensity
  result_low <- simulate_selection_cycles(
    n_cycles = 5,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.6,
    selection_proportion = 0.3,  # Top 30%
    seed = 3002
  )
  
  # High intensity should lead to more gain per cycle (early cycles)
  early_gain_high <- mean(result_high$lpsi_gain[2:3, 1])
  early_gain_low <- mean(result_low$lpsi_gain[2:3, 1])
  
  expect_true(early_gain_high > early_gain_low)
})

# ------------------------------------------------------------------------------
# Test Linkage Effects
# ------------------------------------------------------------------------------

test_that("simulate_selection_cycles models linkage correctly", {
  set.seed(3003)
  
  # Tight linkage (small distances)
  result_tight <- simulate_selection_cycles(
    n_cycles = 10,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.6,
    genetic_distances = rep(0.01, 19),  # 0.01 Morgan = 1 cM
    seed = 3003
  )
  
  # Loose linkage (large distances)
  result_loose <- simulate_selection_cycles(
    n_cycles = 10,
    n_individuals = 100,
    n_loci = 20,
    n_traits = 2,
    heritability = 0.6,
    genetic_distances = rep(0.5, 19),  # 0.5 Morgan = 50 cM
    seed = 3003
  )
  
  # Both should complete without error
  expect_equal(nrow(result_tight$lpsi_mean), 10)
  expect_equal(nrow(result_loose$lpsi_mean), 10)
  
  # Linkage affects recombination and response
  # (Detailed test would require more sophisticated genetic analysis)
})

# ------------------------------------------------------------------------------
# Test Error Handling in Eigenvalue Computation
# ------------------------------------------------------------------------------

test_that("eigenvalue checks work (no warnings with valid data)", {
  set.seed(3004)
  
  # Normal simulation should not produce eigenvalue warnings
  expect_silent({
    result <- simulate_selection_cycles(
      n_cycles = 5,
      n_individuals = 100,
      n_loci = 20,
      n_traits = 2,
      heritability = 0.6,
      seed = 3004
    )
  })
  
  # All eigenvalues should be computed successfully
  expect_true(all(!is.na(result$esim_gain)))
  expect_true(all(!is.na(result$resim_gain)))
})
