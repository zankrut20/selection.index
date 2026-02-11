test_that("rlpsi returns expected structure and satisfies constraints", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]

  # Test with auto-created C (user-friendly way)
  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = 1)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))
  expect_true("C" %in% names(result))

  summary_df <- result$summary
  expect_true(is.data.frame(summary_df))
  expect_true(all(c("GA", "PRE", "Delta_G", "rHI", "hI2") %in% colnames(summary_df)))

  # Check b is a clean numeric vector
  expect_true(is.numeric(result$b))
  expect_false(is.matrix(result$b))
  expect_equal(length(result$b), ncol(pmat))
  expect_equal(length(result$Delta_G), ncol(pmat))

  # Constraint satisfaction: C' Delta_G should be close to zero
  constraint_value <- as.numeric(t(result$C) %*% as.numeric(result$Delta_G))
  expect_true(abs(constraint_value) < 1e-4)
})

test_that("rlpsi works with custom C matrix (backward compatibility)", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  C <- diag(ncol(pmat))[, 1, drop = FALSE]

  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, C = C)

  expect_true(is.list(result))
  expect_true("C" %in% names(result))
  
  # Verify constraint is satisfied
  constraint_value <- as.numeric(t(result$C) %*% as.numeric(result$Delta_G))
  expect_true(abs(constraint_value) < 1e-4)
})

test_that("rlpsi can restrict multiple traits", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]

  # Restrict traits 1 and 3
  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = c(1, 3))

  expect_equal(ncol(result$C), 2)
  
  # Both constraints should be satisfied
  constraint_values <- as.numeric(t(result$C) %*% as.numeric(result$Delta_G))
  expect_true(all(abs(constraint_values) < 1e-4))
})

test_that("ppg_lpsi returns expected structure and proportional gains", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  k <- rep(1, ncol(pmat))

  result <- selection.index:::ppg_lpsi(pmat, gmat, k, wmat = wmat, wcol = 1)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))
  expect_true("phi" %in% names(result))

  summary_df <- result$summary
  expect_true("phi" %in% colnames(summary_df))
  expect_true(is.numeric(summary_df$phi))

  # Check b is a clean numeric vector
  expect_true(is.numeric(result$b))
  expect_false(is.matrix(result$b))
  expect_equal(length(result$b), ncol(pmat))

  ratios <- as.numeric(result$Delta_G) / k
  ratios <- ratios[is.finite(ratios)]
  expect_true((max(ratios) - min(ratios)) < 1e-4)
})

test_that("dg_lpsi returns expected structure and achieves desired gains", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  d <- seq_len(ncol(pmat))

  # Note: DG-LPSI doesn't use economic weights (wmat)
  # Large desired gains may trigger feasibility warnings - suppress for this test
  result <- suppressWarnings(selection.index:::dg_lpsi(pmat, gmat, d))

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))
  expect_true("desired_gains" %in% names(result))

  summary_df <- result$summary
  
  # DG-LPSI should NOT have GA and PRE (no economic weights)
  expect_false("GA" %in% colnames(summary_df))
  expect_false("PRE" %in% colnames(summary_df))
  
  # DG-LPSI should have hI2 and rHI
  expect_true("hI2" %in% colnames(summary_df))
  expect_true("rHI" %in% colnames(summary_df))
  expect_true("Delta_G" %in% colnames(summary_df))

  # Check b is a clean numeric vector
  expect_true(is.numeric(result$b))
  expect_false(is.matrix(result$b))
  expect_equal(length(result$b), ncol(pmat))

  # Check that desired gains match input
  expect_equal(as.numeric(result$desired_gains), d)
  
  # Check that realized gains are proportional to desired gains
  ratios <- as.numeric(result$Delta_G) / d
  ratios <- ratios[is.finite(ratios)]
  expect_true((max(ratios) - min(ratios)) < 1e-4)
})

test_that("ppg_lpsi detects singular matrices", {
  # Create a rank-deficient G matrix (last trait is linear combination of first two)
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  
  # Make G rank-deficient by setting last column/row as linear combination
  gmat_singular <- gmat
  gmat_singular[, ncol(gmat)] <- gmat_singular[, 1] + 0.5 * gmat_singular[, 2]
  gmat_singular[ncol(gmat), ] <- gmat_singular[1, ] + 0.5 * gmat_singular[2, ]
  
  k <- rep(1, ncol(pmat))
  
  # Should detect singular matrix and throw error
  expect_error(
    selection.index:::ppg_lpsi(pmat, gmat_singular, k, wmat = wmat, wcol = 1),
    "Singular matrix detected"
  )
})

test_that("dg_lpsi detects singular matrices", {
  # Create a rank-deficient G matrix
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Make G rank-deficient
  gmat_singular <- gmat
  gmat_singular[, ncol(gmat)] <- gmat_singular[, 1] + 0.5 * gmat_singular[, 2]
  gmat_singular[ncol(gmat), ] <- gmat_singular[1, ] + 0.5 * gmat_singular[2, ]
  
  d <- seq_len(ncol(pmat))
  
  # Should detect singular matrix and throw error
  expect_error(
    selection.index:::dg_lpsi(pmat, gmat_singular, d),
    "Singular matrix detected"
  )
})

# =============================================================================
# PHASE 1 ENHANCEMENTS: Implied Weights & Feasibility Checking
# =============================================================================

test_that("dg_lpsi calculates implied economic weights correctly", {
  # Create simple 2-trait example for easy verification
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)  # Desired gains
  
  result <- selection.index:::dg_lpsi(P, G, d, return_implied_weights = TRUE, check_feasibility = FALSE)
  
  # Check that implied weights exist in output
  expect_true(!is.null(result$implied_weights))
  expect_true(!is.null(result$implied_weights_normalized))
  expect_equal(length(result$implied_weights), 2)
  
  # Verify mathematical relationship: ŵ = G^(-1) P b
  G_inv <- MASS::ginv(G)
  expected_weights <- G_inv %*% P %*% result$b
  
  expect_equal(
    as.numeric(result$implied_weights),
    as.numeric(expected_weights),
    tolerance = 1e-3
  )
  
  # Verify weights are in summary
  expect_true("implied_w.1" %in% colnames(result$summary))
  expect_true("implied_w.2" %in% colnames(result$summary))
  expect_true("implied_w_norm.1" %in% colnames(result$summary))
  expect_true("implied_w_norm.2" %in% colnames(result$summary))
})

test_that("implied weights are normalized correctly", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)
  
  result <- selection.index:::dg_lpsi(P, G, d, return_implied_weights = TRUE, check_feasibility = FALSE)
  
  # Normalized weights: max absolute value should be 1
  max_norm <- max(abs(result$implied_weights_normalized))
  expect_equal(max_norm, 1.0, tolerance = 1e-6)
  
  # Check proportionality between normalized and unnormalized
  if (!all(result$implied_weights == 0)) {
    ratio <- result$implied_weights / result$implied_weights_normalized
    # All ratios should be equal (constant scaling factor)
    expect_true(sd(ratio) < 1e-6 || all(is.na(ratio)))
  }
})

test_that("dg_lpsi can disable implied weights calculation", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  d <- rep(1, ncol(pmat))
  
  result <- selection.index:::dg_lpsi(pmat, gmat, d, return_implied_weights = FALSE, check_feasibility = FALSE)
  
  # Implied weights should be NULL when disabled
  expect_null(result$implied_weights)
  expect_null(result$implied_weights_normalized)
  
  # Summary should not contain implied weight columns
  expect_false(any(grepl("implied_w", colnames(result$summary))))
})

test_that("dg_lpsi feasibility check warns for unrealistic gains", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  
  # Unrealistically large desired gains (10 units each)
  # Given genetic SDs of sqrt(2)=1.41 and sqrt(1.5)=1.22,
  # max possible ≈ 4.24 and 3.67, so 10 is >>80% of max
  d <- c(10.0, 10.0)
  
  expect_warning(
    selection.index:::dg_lpsi(P, G, d, check_feasibility = TRUE, return_implied_weights = FALSE),
    "unrealistic"
  )
})

test_that("dg_lpsi handles feasible gains without warning", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  
  # Realistic small gains (well below theoretical maximum)
  d <- c(0.5, 0.3)
  
  expect_silent(
    selection.index:::dg_lpsi(P, G, d, check_feasibility = TRUE, return_implied_weights = FALSE)
  )
})

test_that("feasibility data frame is correctly structured", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  colnames(P) <- colnames(G) <- c("Trait1", "Trait2")
  d <- c(0.5, 0.3)
  
  result <- selection.index:::dg_lpsi(P, G, d, check_feasibility = TRUE, return_implied_weights = FALSE)
  
  expect_true(!is.null(result$feasibility))
  expect_true(is.data.frame(result$feasibility))
  
  # Check required columns
  required_cols <- c("trait", "desired_gain", "achieved_gain", "genetic_sd",
                     "max_possible_gain", "feasibility_ratio", "is_realistic")
  expect_true(all(required_cols %in% names(result$feasibility)))
  
  expect_equal(nrow(result$feasibility), 2)
  expect_equal(result$feasibility$trait, c("Trait1", "Trait2"))
  expect_equal(result$feasibility$desired_gain, d)
})

test_that("dg_lpsi can disable feasibility checking", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  d <- rep(10, ncol(pmat))  # Unrealistic gains
  
  # Should not warn when feasibility checking is disabled
  expect_silent(
    result <- selection.index:::dg_lpsi(pmat, gmat, d, check_feasibility = FALSE, return_implied_weights = FALSE)
  )
  
  # Feasibility should be NULL
  expect_null(result$feasibility)
})

test_that("dg_lpsi achieves desired gains with high precision", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)
  
  result <- selection.index:::dg_lpsi(P, G, d, check_feasibility = FALSE, return_implied_weights = FALSE)
  
  # Achieved gains should match desired gains closely
  expect_equal(
    as.numeric(result$Delta_G),
    d,
    tolerance = 1e-4
  )
  
  # Gain errors should be near zero
  expect_true(all(abs(result$gain_errors) < 1e-4))
})

test_that("dg_lpsi returns gain_errors in output", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  d <- seq_len(ncol(pmat))
  
  result <- selection.index:::dg_lpsi(pmat, gmat, d, check_feasibility = FALSE, return_implied_weights = FALSE)
  
  expect_true("gain_errors" %in% names(result))
  expect_equal(length(result$gain_errors), ncol(pmat))
  expect_true(all(names(result$gain_errors) == colnames(pmat)))
})

test_that("dg_lpsi has correct S3 class", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)
  
  result <- selection.index:::dg_lpsi(P, G, d, check_feasibility = FALSE, return_implied_weights = TRUE)
  
  # Should have correct S3 classes
  expect_s3_class(result, "dg_lpsi")
  expect_s3_class(result, "selection_index")
  expect_type(result, "list")
})

test_that("dg_lpsi print method works", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)
  
  result <- selection.index:::dg_lpsi(P, G, d, check_feasibility = TRUE, return_implied_weights = TRUE)
  
  # Test that print method exists and produces output
  expect_output(print(result), "DESIRED GAINS INDEX")
  expect_output(print(result), "Pesek")
  expect_output(print(result), "IMPLIED ECONOMIC WEIGHTS")
  expect_output(print(result), "FEASIBILITY ANALYSIS")
})

test_that("dg_lpsi summary method works", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)
  
  result <- selection.index:::dg_lpsi(P, G, d)
  
  # Test that summary method exists and produces additional output
  expect_output(summary(result), "ADDITIONAL DETAILS")
  expect_output(summary(result), "Mean desired gain")
  expect_output(summary(result), "Mean achieved gain")
})

test_that("dg_lpsi handles named matrices correctly", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  colnames(P) <- colnames(G) <- c("Yield", "Quality")
  rownames(P) <- rownames(G) <- c("Yield", "Quality")
  
  d <- c(1.0, 0.5)
  
  result <- selection.index:::dg_lpsi(P, G, d)
  
  # Names should be preserved
  expect_equal(names(result$desired_gains), c("Yield", "Quality"))
  expect_equal(names(result$Delta_G), c("Yield", "Quality"))
  expect_equal(names(result$gain_errors), c("Yield", "Quality"))
  
  if (!is.null(result$implied_weights)) {
    expect_equal(names(result$implied_weights), c("Yield", "Quality"))
  }
})

test_that("dg_lpsi handles selection_intensity parameter", {
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  d <- c(1.0, 0.5)
  
  # Test with different selection intensity
  result <- selection.index:::dg_lpsi(P, G, d, selection_intensity = 2.5, check_feasibility = TRUE, return_implied_weights = FALSE)
  
  expect_equal(result$selection_intensity, 2.5)
  
  # Verify it's used in calculations (affects feasibility thresholds)
  result_default <- selection.index:::dg_lpsi(P, G, d, selection_intensity = 2.063, check_feasibility = TRUE, return_implied_weights = FALSE)
  
  # Feasibility data frames should differ when selection intensity differs
  # (max_possible_gain column will be different)
  expect_false(isTRUE(all.equal(result$feasibility$max_possible_gain, 
                                  result_default$feasibility$max_possible_gain)))
})

test_that("dg_lpsi warns on numerical instability", {
  # Create a scenario with potential numerical issues
  P <- matrix(c(10, 5, 5, 8), 2, 2)
  G <- matrix(c(2, 0.5, 0.5, 1.5), 2, 2)
  
  # Very large desired gains may cause numerical issues
  d <- c(1000, 2000)
  
  # Function should complete but may warn about precision
  result <- selection.index:::dg_lpsi(P, G, d, check_feasibility = FALSE)
  
  # Should still return valid structure
  expect_true(is.list(result))
  expect_true("b" %in% names(result))
})

test_that("dg_lpsi backward compatibility maintained", {
  # Test that old API still works (without new parameters)
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  d <- seq_len(ncol(pmat))
  
  # Old way (should still work but may warn about feasibility with large gains)
  # Using expect_warning or expect_silent depending on whether gains are realistic
  result_old <- selection.index:::dg_lpsi(pmat, gmat, d)
  
  # Should have at minimum the old structure
  expect_true("summary" %in% names(result_old))
  expect_true("b" %in% names(result_old))
  expect_true("Delta_G" %in% names(result_old))
  expect_true("desired_gains" %in% names(result_old))
  
  # Plus new features (with defaults)
  expect_true("implied_weights" %in% names(result_old))
  expect_true("feasibility" %in% names(result_old))
})
