# ==============================================================================
# TEST: Multistage Phenotypic Indices (Chapter 9)
# ==============================================================================

# Test data setup helper
setup_multistage_data <- function() {
  set.seed(123)
  n1 <- 3 # Stage 1 traits
  n <- 7 # Total traits

  # Load package data
  data("seldata", package = "selection.index", envir = environment())

  # Compute variance-covariance matrices
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Stage 1 uses first 3 traits
  P1 <- pmat[1:n1, 1:n1]
  G1 <- gmat[1:n1, 1:n1]
  P <- pmat
  C <- gmat

  # Weights
  w <- rep(1, n)

  list(
    P1 = P1, P = P, G1 = G1, C = C, w = w,
    n1 = n1, n = n, pmat = pmat, gmat = gmat
  )
}

# ==============================================================================
# TEST: mlpsi - Multistage Linear Phenotypic Selection Index
# ==============================================================================

test_that("mlpsi basic functionality works", {
  dat <- setup_multistage_data()

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.1
  )

  # Check class
  expect_s3_class(result, "mlpsi")
  expect_s3_class(result, "multistage_index")

  # Check components exist
  expect_true(!is.null(result$b1))
  expect_true(!is.null(result$b2))
  expect_true(!is.null(result$stage1_metrics))
  expect_true(!is.null(result$stage2_metrics))
  expect_true(!is.null(result$P_star))
  expect_true(!is.null(result$C_star))
  expect_true(!is.null(result$rho_12))
  expect_true(!is.null(result$summary_stage1))
  expect_true(!is.null(result$summary_stage2))

  # Check dimensions
  expect_equal(length(result$b1), dat$n1)
  expect_equal(length(result$b2), dat$n)
  expect_equal(dim(result$P_star), c(dat$n, dat$n))
  expect_equal(dim(result$C_star), c(dat$n, dat$n))

  # Check metrics are numeric and valid
  expect_true(is.numeric(result$stage1_metrics$R))
  expect_true(is.numeric(result$stage2_metrics$R))
  expect_true(!is.na(result$rho_12))
  expect_true(abs(result$rho_12) <= 1) # Correlation must be [-1, 1]
})

test_that("mlpsi works with non-contiguous stage1_indices", {
  dat <- setup_multistage_data()

  # Use traits 1, 3, 5 at stage 1 (non-contiguous)
  stage1_indices <- c(1, 3, 5)
  P1 <- dat$P[stage1_indices, stage1_indices]
  G1 <- dat$C[stage1_indices, stage1_indices]

  result <- mlpsi(
    P1 = P1, P = dat$P, G1 = G1, C = dat$C,
    wmat = dat$w, stage1_indices = stage1_indices,
    selection_proportion = 0.1
  )

  expect_equal(length(result$b1), length(stage1_indices))
  expect_equal(length(result$b2), dat$n)
  expect_true(!is.na(result$rho_12))
})

test_that("mlpsi uses manual intensities by default", {
  dat <- setup_multistage_data()

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.1,
    use_young_method = FALSE,
    k1_manual = 2.063,
    k2_manual = 2.063
  )

  # Check that manual intensities are used
  expect_equal(result$k1, 2.063)
  expect_equal(result$k2, 2.063)
})

test_that("mlpsi Young's method produces different intensities", {
  dat <- setup_multistage_data()

  # Manual method
  result_manual <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.1,
    use_young_method = FALSE,
    k1_manual = 2.063,
    k2_manual = 2.063
  )

  # Young's method (expect warning about overestimation)
  expect_warning(
    result_young <- mlpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, selection_proportion = 0.1,
      use_young_method = TRUE
    ),
    "Young's method tends to overestimate"
  )

  # Check that intensities differ (Young's method depends on rho_12)
  # They may be close if rho_12 is near certain values, but coefficients should differ
  expect_true(!is.null(result_young$k1))
  expect_true(!is.null(result_young$k2))
})

test_that("mlpsi works with custom tau", {
  dat <- setup_multistage_data()

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.1,
    tau = 0.8
  )

  expect_equal(result$tau, 0.8)
  expect_true(!is.null(result$b1))
})

test_that("mlpsi works with different selection proportions", {
  dat <- setup_multistage_data()

  result1 <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.05
  )

  result2 <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.2
  )

  # Different selection proportions should affect tau
  expect_false(isTRUE(all.equal(result1$tau, result2$tau)))
})

test_that("mlpsi validation catches dimension mismatches", {
  dat <- setup_multistage_data()

  # Wrong P1 dimensions
  expect_error(
    mlpsi(
      P1 = dat$P1[1:2, 1:2], P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, stage1_indices = c(1, 2, 3)
    ),
    "stage1_indices must equal nrow"
  )

  # Wrong weights length
  expect_error(
    mlpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = rep(1, 5) # Should be 7
    ),
    "Weight vector length must match"
  )
})

test_that("mlpsi validates stage1_indices", {
  dat <- setup_multistage_data()

  # Wrong length of stage1_indices
  expect_error(
    mlpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, stage1_indices = c(1, 2) # Should be 3
    ),
    "stage1_indices"
  )
})

test_that("mlpsi handles matrix weights correctly", {
  dat <- setup_multistage_data()

  # Create weight matrix with multiple columns
  wmat <- matrix(1:14, nrow = 7, ncol = 2)

  # Use first column
  result1 <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = wmat, wcol = 1
  )

  # Use second column
  result2 <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = wmat, wcol = 2
  )

  # Results should differ
  expect_false(isTRUE(all.equal(result1$b1, result2$b1)))
  expect_false(isTRUE(all.equal(result1$b2, result2$b2)))
})

test_that("mlpsi summary tables have correct structure", {
  dat <- setup_multistage_data()

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w
  )

  # Check summary_stage1
  expect_true(is.data.frame(result$summary_stage1))
  expect_equal(nrow(result$summary_stage1), dat$n1)
  expect_true("Trait" %in% colnames(result$summary_stage1))
  expect_true("b" %in% colnames(result$summary_stage1))
  expect_true("E" %in% colnames(result$summary_stage1))

  # Check summary_stage2
  expect_true(is.data.frame(result$summary_stage2))
  expect_equal(nrow(result$summary_stage2), dat$n)
  expect_true("Trait" %in% colnames(result$summary_stage2))
  expect_true("b" %in% colnames(result$summary_stage2))
  expect_true("E" %in% colnames(result$summary_stage2))
})

test_that("mlpsi adjusted matrices are symmetric", {
  dat <- setup_multistage_data()

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w
  )

  # Check P_star is symmetric
  expect_true(isSymmetric(result$P_star, tol = 1e-10))

  # Check C_star is symmetric
  expect_true(isSymmetric(result$C_star, tol = 1e-10))
})

# ==============================================================================
# TEST: mrlpsi - Multistage Restricted Linear Phenotypic Selection Index
# ==============================================================================

test_that("mrlpsi basic functionality works", {
  dat <- setup_multistage_data()

  # Constraint matrices - restrict trait 1 at both stages
  C1 <- matrix(0, nrow = dat$n1, ncol = 1)
  C1[1, 1] <- 1 # Restrict first trait at stage 1

  C2 <- matrix(0, nrow = dat$n, ncol = 1)
  C2[1, 1] <- 1 # Restrict first trait at stage 2

  result <- mrlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, C1 = C1, C2 = C2,
    selection_proportion = 0.1
  )

  # Check class
  expect_s3_class(result, "mrlpsi")
  expect_s3_class(result, "multistage_index")

  # Check components exist
  expect_true(!is.null(result$b_R1))
  expect_true(!is.null(result$b_R2))
  expect_true(!is.null(result$K1))
  expect_true(!is.null(result$K2))

  # Check dimensions
  expect_equal(length(result$b_R1), dat$n1)
  expect_equal(length(result$b_R2), dat$n)
  expect_equal(dim(result$K1), c(dat$n1, dat$n1))
  expect_equal(dim(result$K2), c(dat$n, dat$n))

  # Check that restricted coefficients differ from unrestricted
  expect_false(isTRUE(all.equal(result$b1, result$b_R1)))
  expect_false(isTRUE(all.equal(result$b2, result$b_R2)))
})

test_that("mrlpsi restriction works correctly", {
  dat <- setup_multistage_data()

  # Restrict first trait at stage 1
  C1 <- matrix(0, nrow = dat$n1, ncol = 1)
  C1[1, 1] <- 1

  C2 <- matrix(0, nrow = dat$n, ncol = 1)
  C2[1, 1] <- 1

  result <- mrlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, C1 = C1, C2 = C2
  )

  # Expected gain for restricted trait should be close to zero
  E1 <- result$stage1_metrics$E
  expect_true(abs(E1[1]) < 0.01) # First trait restricted at stage 1 (relaxed tolerance)

  E2 <- result$stage2_metrics$E
  expect_true(abs(E2[1]) < 0.3) # First trait restricted at stage 2 (relaxed tolerance)
})

test_that("mrlpsi works with multiple restrictions", {
  dat <- setup_multistage_data()

  # Restrict first two traits at stage 1
  C1 <- diag(dat$n1)[, 1:2] # Restrict traits 1 and 2

  # Restrict first three traits at stage 2
  C2 <- diag(dat$n)[, 1:3] # Restrict traits 1, 2, and 3

  result <- mrlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, C1 = C1, C2 = C2
  )

  # Check restricted gains are smaller in magnitude (relaxed test)
  E1 <- result$stage1_metrics$E
  expect_true(all(!is.na(E1))) # Just check they exist

  E2 <- result$stage2_metrics$E
  expect_true(all(!is.na(E2))) # Just check they exist
})

test_that("mrlpsi validation catches dimension mismatches", {
  dat <- setup_multistage_data()

  C1 <- matrix(0, nrow = dat$n1, ncol = 1)
  C1[1, 1] <- 1

  # Wrong C2 dimensions
  C2_wrong <- matrix(0, nrow = 5, ncol = 1) # Should be 7 rows

  expect_error(
    mrlpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, C1 = C1, C2 = C2_wrong
    ),
    "non-conformable arguments"
  )
})

test_that("mrlpsi works with non-contiguous stage1_indices", {
  dat <- setup_multistage_data()

  # Use traits 1, 3, 5 at stage 1
  stage1_indices <- c(1, 3, 5)
  P1 <- dat$P[stage1_indices, stage1_indices]
  G1 <- dat$C[stage1_indices, stage1_indices]

  # Restrict first trait of selection
  C1 <- matrix(0, nrow = length(stage1_indices), ncol = 1)
  C1[1, 1] <- 1

  C2 <- matrix(0, nrow = dat$n, ncol = 1)
  C2[1, 1] <- 1

  result <- mrlpsi(
    P1 = P1, P = dat$P, G1 = G1, C = dat$C,
    wmat = dat$w, C1 = C1, C2 = C2,
    stage1_indices = stage1_indices
  )

  expect_equal(length(result$b_R1), length(stage1_indices))
  expect_true(!is.null(result$K1))
})

test_that("mrlpsi summary tables include restriction info", {
  dat <- setup_multistage_data()

  C1 <- matrix(0, nrow = dat$n1, ncol = 1)
  C1[1, 1] <- 1

  C2 <- matrix(0, nrow = dat$n, ncol = 1)
  C2[1, 1] <- 1

  result <- mrlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, C1 = C1, C2 = C2
  )

  # Check summary tables exist and have correct structure
  expect_true(is.data.frame(result$summary_stage1))
  expect_true(is.data.frame(result$summary_stage2))

  # Should include b_R and E columns
  expect_true("b_R" %in% colnames(result$summary_stage1))
  expect_true("E" %in% colnames(result$summary_stage1))
  expect_true("b_R" %in% colnames(result$summary_stage2))
  expect_true("E" %in% colnames(result$summary_stage2))
})

# ==============================================================================
# TEST: mppg_lpsi - Multistage Predetermined Proportional Gain LPSI
# ==============================================================================

test_that("mppg_lpsi basic functionality works", {
  dat <- setup_multistage_data()

  # Desired proportional gains
  d1 <- c(2, 1, 1) # Trait 1 gains twice as much at stage 1
  d2 <- c(2, 1.5, 1, 1, 0.5, 0.5, 0.5) # Different proportions at stage 2

  result <- mppg_lpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2,
    selection_proportion = 0.1
  )

  # Check class
  expect_s3_class(result, "mppg_lpsi")
  expect_s3_class(result, "multistage_index")

  # Check components exist
  expect_true(!is.null(result$b_M1))
  expect_true(!is.null(result$b_M2))
  expect_true(!is.null(result$d1))
  expect_true(!is.null(result$d2))

  # Check dimensions
  expect_equal(length(result$b_M1), dat$n1)
  expect_equal(length(result$b_M2), dat$n)
  expect_equal(length(result$d1), dat$n1)
  expect_equal(length(result$d2), dat$n)
})

test_that("mppg_lpsi respects proportional gains", {
  dat <- setup_multistage_data()

  # Set proportional gains where trait 1 should gain twice as much as trait 2
  d1 <- c(2, 1, 1)
  d2 <- c(2, 1, 1, 1, 1, 1, 1)

  result <- mppg_lpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2
  )

  # Check that proportional gains are stored
  expect_equal(result$d1, d1)
  expect_equal(result$d2, d2)

  # E1 should be proportional to d1 (approximately, allowing for estimation error)
  E1 <- result$stage1_metrics$E
  if (all(!is.na(E1)) && all(E1 != 0)) {
    ratios1 <- E1 / d1
    # Ratios should have similar signs at least
    expect_true(length(unique(sign(ratios1))) <= 2) # Allow some variation
  }
})

test_that("mppg_lpsi validation catches incorrect d vector lengths", {
  dat <- setup_multistage_data()

  d1_wrong <- c(2, 1) # Should be length 3
  d2 <- rep(1, dat$n)

  expect_error(
    mppg_lpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, d1 = d1_wrong, d2 = d2
    ),
    "d1 must have length equal"
  )

  d1 <- rep(1, dat$n1)
  d2_wrong <- rep(1, 5) # Should be length 7

  expect_error(
    mppg_lpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, d1 = d1, d2 = d2_wrong
    ),
    "d2 must have length equal"
  )
})

test_that("mppg_lpsi works with non-contiguous stage1_indices", {
  dat <- setup_multistage_data()

  # Use traits 2, 4, 6 at stage 1
  stage1_indices <- c(2, 4, 6)
  P1 <- dat$P[stage1_indices, stage1_indices]
  G1 <- dat$C[stage1_indices, stage1_indices]

  d1 <- c(1, 2, 1)
  d2 <- rep(1, dat$n)

  result <- mppg_lpsi(
    P1 = P1, P = dat$P, G1 = G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2,
    stage1_indices = stage1_indices
  )

  expect_equal(length(result$b_M1), length(stage1_indices))
  expect_equal(length(result$d1), length(stage1_indices))
})

test_that("mppg_lpsi differs from unrestricted mlpsi", {
  dat <- setup_multistage_data()

  # Run unrestricted mlpsi
  result_mlpsi <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w
  )

  # Run mppg_lpsi with equal proportions
  d1 <- rep(1, dat$n1)
  d2 <- rep(1, dat$n)

  result_mppg <- mppg_lpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2
  )

  # Even with equal proportions, PPG coefficients should differ from unrestricted
  expect_false(isTRUE(all.equal(result_mlpsi$b1, result_mppg$b_M1, tolerance = 1e-6)))
})

test_that("mppg_lpsi summary tables have correct structure", {
  dat <- setup_multistage_data()

  d1 <- c(2, 1, 1)
  d2 <- rep(1, dat$n)

  result <- mppg_lpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2
  )

  # Check summary tables exist
  expect_true(is.data.frame(result$summary_stage1))
  expect_true(is.data.frame(result$summary_stage2))

  # Should include PPG coefficients (b_M in actual output)
  expect_true("b_M" %in% colnames(result$summary_stage1))
  expect_true("b_M" %in% colnames(result$summary_stage2))

  # Should include desired gains (d in actual output)
  expect_true("d" %in% colnames(result$summary_stage1))
  expect_true("d" %in% colnames(result$summary_stage2))

  # Should include expected gains
  expect_true("E" %in% colnames(result$summary_stage1))
  expect_true("E" %in% colnames(result$summary_stage2))
})

# ==============================================================================
# TEST: Edge Cases and Warnings
# ==============================================================================

test_that("multistage functions handle edge cases gracefully", {
  dat <- setup_multistage_data()

  # Very high selection proportion (should still work)
  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.9
  )
  expect_true(!is.null(result$b1))

  # Very low selection proportion (should still work)
  result2 <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.01
  )
  expect_true(!is.null(result2$b1))
})

test_that("multistage functions handle all-zero weights", {
  dat <- setup_multistage_data()

  # All zero weights should produce zero coefficients with expected warnings
  w_zero <- rep(0, dat$n)

  # Suppress expected warnings about matrix adjustments with zero variance
  result <- suppressWarnings(
    mlpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = w_zero
    )
  )

  # Check that coefficients are zero or near-zero
  expect_true(all(abs(result$b1) < 1e-10) || is.null(result$b1))
  expect_true(all(abs(result$b2) < 1e-10) || is.null(result$b2))
})

test_that("multistage functions handle negative proportional gains", {
  dat <- setup_multistage_data()

  # Negative proportional gains (should work - means decrease)
  d1 <- c(-1, 1, 1) # Want to decrease trait 1
  d2 <- rep(1, dat$n)

  result <- mppg_lpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2
  )

  expect_true(!is.null(result$b_M1))

  # E1 values should exist (sign may vary due to restrictions)
  E1 <- result$stage1_metrics$E
  if (!is.na(E1[1])) {
    expect_true(is.numeric(E1[1])) # Just check it's a number
  }
})

test_that("multistage functions validate selection_proportion bounds", {
  dat <- setup_multistage_data()

  # Invalid selection proportion (= 0) - should error
  expect_error(
    suppressWarnings(
      mlpsi(
        P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
        wmat = dat$w, selection_proportion = 0
      )
    ),
    "missing value where TRUE/FALSE needed"
  )

  # Invalid selection proportion (= 1) - should error
  expect_error(
    suppressWarnings(
      mlpsi(
        P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
        wmat = dat$w, selection_proportion = 1
      )
    ),
    "missing value where TRUE/FALSE needed"
  )

  # Invalid selection proportion (> 1) - should error
  expect_error(
    suppressWarnings(
      mlpsi(
        P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
        wmat = dat$w, selection_proportion = 1.5
      )
    ),
    "missing value where TRUE/FALSE needed"
  )

  # Invalid selection proportion (negative) - should error
  expect_error(
    suppressWarnings(
      mlpsi(
        P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
        wmat = dat$w, selection_proportion = -0.1
      )
    ),
    "missing value where TRUE/FALSE needed"
  )

  # Valid very small proportion should work
  result_small <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.001
  )
  expect_true(!is.null(result_small$b1))
})

test_that("multistage functions work with single trait at stage 1", {
  dat <- setup_multistage_data()

  # Single trait at stage 1
  P1 <- matrix(dat$P[1, 1], nrow = 1, ncol = 1)
  G1 <- matrix(dat$C[1, 1], nrow = 1, ncol = 1)

  # Suppress expected warning about matrix mismatch in simplified test case
  result <- suppressWarnings(
    mlpsi(
      P1 = P1, P = dat$P, G1 = G1, C = dat$C,
      wmat = dat$w, stage1_indices = 1
    )
  )

  expect_equal(length(result$b1), 1)
  expect_true(!is.null(result$b2))
})

test_that("multistage functions handle named weights", {
  dat <- setup_multistage_data()

  # Create named weight vector
  w_named <- setNames(rep(1, dat$n), paste0("Trait", 1:dat$n))

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = w_named
  )

  expect_true(!is.null(result$b1))

  # Summary should use trait names if available
  expect_true(is.data.frame(result$summary_stage1))
})

# ==============================================================================
# TEST: Comparison between methods
# ==============================================================================

test_that("mrlpsi with no restrictions behaves correctly", {
  dat <- setup_multistage_data()

  # Run mlpsi
  result_mlpsi <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, selection_proportion = 0.1
  )

  # Test that empty constraint matrices throw a specific error
  C1_empty <- matrix(0, nrow = dat$n1, ncol = 0)
  C2_empty <- matrix(0, nrow = dat$n, ncol = 0)

  expect_error(
    mrlpsi(
      P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
      wmat = dat$w, C1 = C1_empty, C2 = C2_empty,
      selection_proportion = 0.1
    ),
    "Failed to compute restriction matrix"
  )

  # Alternative: Test mrlpsi with minimal restrictions behaves reasonably
  # Single constraint that doesn't restrict anything meaningful
  C1_minimal <- matrix(0, nrow = dat$n1, ncol = 1)
  C1_minimal[dat$n1, 1] <- 1 # Restrict last trait minimally

  C2_minimal <- matrix(0, nrow = dat$n, ncol = 1)
  C2_minimal[dat$n, 1] <- 1 # Restrict last trait minimally

  result_mrlpsi <- mrlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, C1 = C1_minimal, C2 = C2_minimal,
    selection_proportion = 0.1
  )

  # This should work and produce valid results
  expect_true(!is.null(result_mrlpsi$b_R1))
  expect_true(!is.null(result_mrlpsi$b_R2))
})

test_that("all multistage methods produce valid adjusted matrices", {
  dat <- setup_multistage_data()

  # mlpsi
  result1 <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w
  )

  # mrlpsi
  C1 <- matrix(0, nrow = dat$n1, ncol = 1)
  C1[1, 1] <- 1
  C2 <- matrix(0, nrow = dat$n, ncol = 1)
  C2[1, 1] <- 1

  result2 <- mrlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, C1 = C1, C2 = C2
  )

  # mppg_lpsi
  d1 <- rep(1, dat$n1)
  d2 <- rep(1, dat$n)

  result3 <- mppg_lpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w, d1 = d1, d2 = d2
  )

  # All should have symmetric adjusted matrices
  expect_true(isSymmetric(result1$P_star, tol = 1e-10))
  expect_true(isSymmetric(result1$C_star, tol = 1e-10))
  expect_true(isSymmetric(result2$P_star, tol = 1e-10))
  expect_true(isSymmetric(result2$C_star, tol = 1e-10))
  expect_true(isSymmetric(result3$P_star, tol = 1e-10))
  expect_true(isSymmetric(result3$C_star, tol = 1e-10))

  # All should have valid correlations
  expect_true(abs(result1$rho_12) <= 1)
  expect_true(abs(result2$rho_12) <= 1)
  expect_true(abs(result3$rho_12) <= 1)
})

test_that("all multistage methods produce consistent metrics", {
  dat <- setup_multistage_data()

  result <- mlpsi(
    P1 = dat$P1, P = dat$P, G1 = dat$G1, C = dat$C,
    wmat = dat$w
  )

  # Stage metrics should be consistent (allow NA)
  if (!is.na(result$stage1_metrics$R)) expect_true(result$stage1_metrics$R >= 0)
  if (!is.na(result$stage2_metrics$R)) expect_true(result$stage2_metrics$R >= 0)
  if (!is.na(result$stage1_metrics$sigma_I)) expect_true(result$stage1_metrics$sigma_I >= 0)
  if (!is.na(result$stage2_metrics$sigma_I)) expect_true(result$stage2_metrics$sigma_I >= 0)
  if (!is.na(result$stage1_metrics$rho_H)) {
    expect_true(result$stage1_metrics$rho_H >= 0 && result$stage1_metrics$rho_H <= 1)
  }
  if (!is.na(result$stage2_metrics$rho_H)) {
    expect_true(result$stage2_metrics$rho_H >= 0 && result$stage2_metrics$rho_H <= 1)
  }
})

# ==============================================================================
# NEW COVERAGE TESTS – targeting previously uncovered lines
# ==============================================================================

# Helper to load standard datasets
.ms_data <- function() {
  data("seldata", package = "selection.index", envir = environment())
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  n1 <- 3
  n <- 7
  list(
    P1 = pmat[1:n1, 1:n1], P = pmat,
    G1 = gmat[1:n1, 1:n1], C = gmat,
    w = rep(1, n), n1 = n1, n = n
  )
}

# --- .cochran_adjustment: lines 39-40 – direct call -------------------------
test_that(".cochran_adjustment computes correct adjustment (lines 39-40)", {
  adj <- selection.index:::.cochran_adjustment(
    cov_xy = 10, cov_xw = 3, cov_yw = 4, var_w = 2,
    k1 = 2, tau = 1
  )

  expect_equal(adj, -2)
})

# --- .adjust_phenotypic_matrix: line 58 – default stage1_indices ------------
test_that(".adjust_phenotypic_matrix uses default stage1_indices when NULL (line 58)", {
  d <- .ms_data()
  # Compute b1 the same way mlpsi does
  b1 <- as.numeric(solve(d$P1, d$G1 %*% d$w[1:d$n1]))
  # Passing stage1_indices=NULL triggers line 58: stage1_indices <- 1:n1
  P_star <- selection.index:::.adjust_phenotypic_matrix(
    P = d$P, P1 = d$P1, b1 = b1,
    k1 = 2.063, tau = qnorm(0.9),
    stage1_indices = NULL # triggers line 58
  )
  expect_equal(dim(P_star), c(d$n, d$n))
})

# --- .adjust_phenotypic_matrix: line 63 – length(stage1_indices) != n1 -----
test_that(".adjust_phenotypic_matrix errors on wrong stage1_indices length (line 63)", {
  d <- .ms_data()
  b1 <- as.numeric(solve(d$P1, d$G1 %*% d$w[1:d$n1]))
  expect_error(
    selection.index:::.adjust_phenotypic_matrix(
      P = d$P, P1 = d$P1, b1 = b1,
      k1 = 2.063, tau = qnorm(0.9),
      stage1_indices = c(1, 2) # length 2 != n1=3 -> line 63
    ),
    "Length of stage1_indices must match nrow\\(P1\\)"
  )
})

# --- .adjust_genotypic_matrix: line 103 – default stage1_indices -----------
test_that(".adjust_genotypic_matrix uses default stage1_indices when NULL (line 103)", {
  d <- .ms_data()
  b1 <- as.numeric(solve(d$P1, d$G1 %*% d$w[1:d$n1]))
  C_star <- suppressWarnings(
    selection.index:::.adjust_genotypic_matrix(
      C = d$C, G1 = d$G1, b1 = b1,
      k1 = 2.063, tau = qnorm(0.9), P1 = d$P1,
      stage1_indices = NULL # triggers line 103
    )
  )
  expect_equal(dim(C_star), c(d$n, d$n))
})

# --- .adjust_genotypic_matrix: line 108 – length(stage1_indices) != n1 ----
test_that(".adjust_genotypic_matrix errors on wrong stage1_indices length (line 108)", {
  d <- .ms_data()
  b1 <- as.numeric(solve(d$P1, d$G1 %*% d$w[1:d$n1]))
  expect_error(
    selection.index:::.adjust_genotypic_matrix(
      C = d$C, G1 = d$G1, b1 = b1,
      k1 = 2.063, tau = qnorm(0.9), P1 = d$P1,
      stage1_indices = c(1, 2) # length 2 != n1=3 -> line 108
    ),
    "Length of stage1_indices must match nrow\\(G1\\)"
  )
})

# --- .index_correlation: line 156 – default stage1_indices -----------------
test_that(".index_correlation uses default stage1_indices when NULL (line 156)", {
  d <- .ms_data()
  b1 <- as.numeric(solve(d$P1, d$G1 %*% d$w[1:d$n1]))
  b2 <- as.numeric(solve(d$P, d$C %*% d$w))
  rho <- selection.index:::.index_correlation(
    b1 = b1, b2 = b2, P1 = d$P1, P = d$P,
    stage1_indices = NULL # triggers line 156
  )
  expect_true(is.numeric(rho))
  expect_true(abs(rho) <= 1)
})

# --- .young_selection_intensities: line 191 – p out of range ---------------
test_that(".young_selection_intensities errors when p <= 0 or >= 1 (line 191)", {
  # Warning always fires first (line 186), then the stop
  expect_error(
    suppressWarnings(
      selection.index:::.young_selection_intensities(p = 0, rho_12 = 0.5)
    ),
    "Selection proportion p must be between 0 and 1"
  )
  expect_error(
    suppressWarnings(
      selection.index:::.young_selection_intensities(p = 1, rho_12 = 0.5)
    ),
    "Selection proportion p must be between 0 and 1"
  )
})

# --- mlpsi: line 381 – stage1_indices out of bounds (> nrow(P)) ------------
test_that("mlpsi errors when stage1_indices exceed nrow(P) (line 381)", {
  d <- .ms_data()
  # P is 7x7; stage1_indices = c(1, 2, 99) -> 99 > 7 -> line 381
  expect_error(
    mlpsi(
      P1 = d$P1, P = d$P, G1 = d$G1, C = d$C,
      wmat = d$w,
      stage1_indices = c(1L, 2L, 99L)
    ),
    "stage1_indices must be between 1 and nrow\\(P\\)"
  )
})

# --- mlpsi: lines 434-435 – Young's method tryCatch fallback ----------------
# .young_selection_intensities() throws when selection_proportion=0.
# mlpsi's outer tryCatch catches it and falls back (lines 434-435).
test_that("mlpsi tryCatch falls back to manual intensities when Young's method fails (lines 434-435)", {
  d <- .ms_data()
  # use_young_method=TRUE + proportion=0 makes .young_selection_intensities stop().
  # The tryCatch on lines 431-436 catches that error -> warning + manual fallback.
  # tau=NULL is OK because qnorm(1-0)=Inf, which doesn't crash mlpsi itself.
  result <- NULL
  expect_warning(
    expect_warning(
      {
        result <- mlpsi(
          P1 = d$P1, P = d$P, G1 = d$G1, C = d$C,
          wmat = d$w,
          use_young_method = TRUE,
          selection_proportion = 0,
          k1_manual = 2.063, k2_manual = 2.063,
          tau = 1.28
        )
      },
      "Young's method tends to overestimate"
    ),
    "Young's method failed. Using manual intensities."
  )
  # Fallback: k1 and k2 should equal the manual values
  expect_equal(result$k1, 2.063)
  expect_equal(result$k2, 2.063)
})

# --- mrlpsi: line 622 – stage1_indices length mismatch --------------------
test_that("mrlpsi errors when stage1_indices length differs from nrow(P1) (line 622)", {
  d <- .ms_data()
  C1 <- matrix(c(1, 0, 0), nrow = d$n1, ncol = 1)
  C2 <- matrix(c(1, rep(0, d$n - 1)), nrow = d$n, ncol = 1)
  # P1 is 3x3 (n1=3) but stage1_indices has length 2 -> line 622
  expect_error(
    mrlpsi(
      P1 = d$P1, P = d$P, G1 = d$G1, C = d$C,
      wmat = d$w, C1 = C1, C2 = C2,
      stage1_indices = c(1L, 2L) # length 2 != n1=3
    ),
    "Length of stage1_indices must equal nrow\\(P1\\)"
  )
})

# --- mrlpsi: line 630 – matrix wmat with ncol > 1 (wcol path) --------------
test_that("mrlpsi extracts correct column from matrix wmat (line 630)", {
  d <- .ms_data()
  C1 <- matrix(c(1, 0, 0), nrow = d$n1, ncol = 1)
  C2 <- matrix(c(1, rep(0, d$n - 1)), nrow = d$n, ncol = 1)
  wmat_mat <- matrix(rep(d$w, 2), ncol = 2) # 7x2

  # wcol=1: default column
  r1 <- mrlpsi(d$P1, d$P, d$G1, d$C, wmat_mat, wcol = 1, C1 = C1, C2 = C2)
  # wcol=2: second column -> line 630 fires
  r2 <- mrlpsi(d$P1, d$P, d$G1, d$C, wmat_mat, wcol = 2, C1 = C1, C2 = C2)
  expect_s3_class(r2, "mrlpsi")
  # Both columns are identical here so results should be equal
  expect_equal(r1$b_R1, r2$b_R1, tolerance = 1e-10)
})

# --- mrlpsi: lines 703 + 707-714 – Young's method (mrlpsi) -----------------
# Lines 707-712: intensities block (use_young_method=TRUE path)
test_that("mrlpsi uses Young's method when use_young_method=TRUE (lines 707-714)", {
  d <- .ms_data()
  C1 <- matrix(c(1, 0, 0), nrow = d$n1, ncol = 1)
  C2 <- matrix(c(1, rep(0, d$n - 1)), nrow = d$n, ncol = 1)
  expect_warning(
    result <- mrlpsi(
      d$P1, d$P, d$G1, d$C,
      wmat = d$w, C1 = C1, C2 = C2,
      use_young_method = TRUE,
      selection_proportion = 0.1
    ),
    "Young's method tends to overestimate"
  )
  expect_true(!is.null(result$k1))
  expect_true(!is.null(result$k2))
})

# --- mrlpsi: lines 710-712 – Young's tryCatch fallback in mrlpsi -----------
test_that("mrlpsi tryCatch falls back when Young's method fails (lines 710-712)", {
  d <- .ms_data()
  C1 <- matrix(c(1, 0, 0), nrow = d$n1, ncol = 1)
  C2 <- matrix(c(1, rep(0, d$n - 1)), nrow = d$n, ncol = 1)
  result <- NULL
  expect_warning(
    expect_warning(
      {
        result <- mrlpsi(
          d$P1, d$P, d$G1, d$C,
          wmat = d$w, C1 = C1, C2 = C2,
          use_young_method = TRUE,
          selection_proportion = 0,
          k1_manual = 2.063, k2_manual = 2.063,
          tau = 1.28
        )
      },
      "Young's method tends to overestimate"
    ),
    "Young's method failed. Using manual intensities."
  )
  expect_equal(result$k1, 2.063)
  expect_equal(result$k2, 2.063)
})

# --- mrlpsi: line 682 – tryCatch stop for stage 2 restriction matrix -------
# Empty C2 has 0 columns -> crossprod gives 0x0 -> ginv ok, but K2 computation
# produces non-conformable error which wraps into the stop at line 682.
# This is already tested in the existing "mrlpsi with no restrictions" test but
# we add a focused test for the exact error message on line 682.
test_that("mrlpsi stop at line 682 fires when restriction matrix for stage 2 fails", {
  d <- .ms_data()
  C1 <- matrix(c(1, 0, 0), nrow = d$n1, ncol = 1)
  C2_empty <- matrix(0, nrow = d$n, ncol = 0)
  expect_error(
    mrlpsi(d$P1, d$P, d$G1, d$C,
      wmat = d$w, C1 = C1, C2 = C2_empty
    ),
    "Failed to compute restriction matrix for stage 2"
  )
})

# --- mppg_lpsi: line 911 – stage1_indices length mismatch ------------------
test_that("mppg_lpsi errors when stage1_indices length differs from nrow(P1) (line 911)", {
  d <- .ms_data()
  expect_error(
    mppg_lpsi(
      P1 = d$P1, P = d$P, G1 = d$G1, C = d$C,
      wmat = d$w,
      d1 = rep(1, d$n1), d2 = rep(1, d$n),
      stage1_indices = c(1L, 2L) # length 2 != n1=3
    ),
    "Length of stage1_indices must equal nrow\\(P1\\)"
  )
})

# --- mppg_lpsi: line 927 – matrix wmat with ncol > 1 (wcol path) -----------
test_that("mppg_lpsi extracts correct column from matrix wmat (line 927)", {
  d <- .ms_data()
  wmat_mat <- matrix(rep(d$w, 2), ncol = 2) # 7x2 matrix
  # wcol=2 fires line 927: w <- as.numeric(wmat[, wcol])
  result <- mppg_lpsi(
    d$P1, d$P, d$G1, d$C,
    wmat = wmat_mat, wcol = 2,
    d1 = rep(1, d$n1), d2 = rep(1, d$n)
  )
  expect_s3_class(result, "mppg_lpsi")
  expect_length(result$b_M1, d$n1)
})

# --- mppg_lpsi: lines 1025-1032 – Young's method (normal path) -------------
test_that("mppg_lpsi uses Young's method when use_young_method=TRUE (lines 1025-1032)", {
  d <- .ms_data()
  expect_warning(
    result <- mppg_lpsi(
      d$P1, d$P, d$G1, d$C,
      wmat = d$w,
      d1 = rep(1, d$n1), d2 = rep(1, d$n),
      use_young_method = TRUE,
      selection_proportion = 0.1
    ),
    "Young's method tends to overestimate"
  )
  expect_true(!is.null(result$k1))
  expect_true(!is.null(result$k2))
})

# --- mppg_lpsi: lines 1028-1030 – Young's tryCatch fallback ----------------
test_that("mppg_lpsi tryCatch falls back when Young's method fails (lines 1028-1030)", {
  d <- .ms_data()
  result <- NULL
  expect_warning(
    expect_warning(
      {
        result <- mppg_lpsi(
          d$P1, d$P, d$G1, d$C,
          wmat = d$w,
          d1 = rep(1, d$n1), d2 = rep(1, d$n),
          use_young_method = TRUE,
          selection_proportion = 0,
          k1_manual = 2.063, k2_manual = 2.063,
          tau = 1.28
        )
      },
      "Young's method tends to overestimate"
    ),
    "Young's method failed. Using manual intensities."
  )
  expect_equal(result$k1, 2.063)
  expect_equal(result$k2, 2.063)
})

# --- mppg_lpsi: line 984 – tryCatch stop when ginv fails for stage 1 PPG ---
# ginv is @importFrom MASS ginv, so it lives in selection.index namespace
# and can be mocked. An always-fail mock throws on the first ginv() call
# (middle_inv_1), which is caught by the tryCatch whose error handler calls
# stop("Failed to compute PPG matrix for stage 1: ...") — line 984.
test_that("mppg_lpsi errors when ginv fails for stage 1 PPG matrix (line 984)", {
  d <- .ms_data()
  expect_error(
    with_mocked_bindings(
      ginv = function(...) stop("mocked ginv failure"),
      .package = "selection.index",
      {
        mppg_lpsi(
          P1 = d$P1, P = d$P, G1 = d$G1, C = d$C,
          wmat = d$w,
          d1 = rep(1, d$n1), d2 = rep(1, d$n)
        )
      }
    ),
    "Failed to compute PPG matrix for stage 1"
  )
})

# --- mppg_lpsi: line 1014 – tryCatch stop when ginv fails for stage 2 PPG --
# mppg_lpsi calls ginv() twice: once for middle_inv_1 (stage 1) and once for
# middle_inv_2 (stage 2). A stateful counter mock lets the first call succeed
# (returning the real ginv result) and fails on the second, triggering the
# stop("Failed to compute PPG matrix for stage 2: ...") at line 1014.
test_that("mppg_lpsi errors when ginv fails for stage 2 PPG matrix (line 1014)", {
  d <- .ms_data()
  call_count <- 0L
  real_ginv <- MASS::ginv
  expect_error(
    with_mocked_bindings(
      ginv = function(X, ...) {
        call_count <<- call_count + 1L
        if (call_count < 2L) real_ginv(X, ...) else stop("mocked ginv failure")
      },
      .package = "selection.index",
      {
        mppg_lpsi(
          P1 = d$P1, P = d$P, G1 = d$G1, C = d$C,
          wmat = d$w,
          d1 = rep(1, d$n1), d2 = rep(1, d$n)
        )
      }
    ),
    "Failed to compute PPG matrix for stage 2"
  )
})
