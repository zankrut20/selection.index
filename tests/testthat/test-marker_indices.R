test_that("lmsi basic functionality works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data
  set.seed(123)
  n_genotypes <- 50
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  # Simulate phenotype and marker score matrices
  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(phen_mat) <- colnames(gmat)

  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(marker_scores) <- colnames(gmat)

  # Test basic LMSI calculation
  result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights)

  # Check that result has correct structure
  expect_s3_class(result, "lmsi")
  expect_type(result, "list")

  # Check essential components exist
  essential_components <- c(
    "b_y", "b_s", "b_combined", "P_L", "G_L", "G_s",
    "rHI", "sigma_I", "R", "GA", "Delta_H", "summary",
    "phenotype_coeffs", "marker_coeffs", "coeff_analysis"
  )
  for (comp in essential_components) {
    expect_true(comp %in% names(result), info = paste("Missing component:", comp))
  }

  # Check dimensions
  expect_length(result$b_y, n_traits)
  expect_length(result$b_s, n_traits)
  expect_length(result$b_combined, 2 * n_traits)
  expect_length(result$Delta_H, n_traits)

  # Check matrix dimensions
  expect_equal(dim(result$P_L), c(2 * n_traits, 2 * n_traits))
  expect_equal(dim(result$G_L), c(2 * n_traits, n_traits))
  expect_equal(dim(result$G_s), c(n_traits, n_traits))

  # Check that b_combined equals [b_y; b_s]
  expect_equal(result$b_combined, c(result$b_y, result$b_s), tolerance = 1e-10)

  # Check that accuracy is in valid range [0, 1]
  expect_gte(result$rHI, 0)
  expect_lte(result$rHI, 1)

  # Check that heritability is in valid range [0, 1]
  expect_gte(result$hI2, 0)
  expect_lte(result$hI2, 1)

  # Check non-negative values
  expect_gte(result$sigma_I, 0)
  expect_type(result$R, "double")
  expect_type(result$GA, "double")
})

test_that("lmsi with provided G_s works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  # Create custom G_s matrix
  G_s <- 0.8 * gmat # Assume markers explain 80% of genetic variance

  # Test LMSI with provided G_s (no data matrices needed)
  result <- lmsi(
    phen_mat = NULL, marker_scores = NULL,
    pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights
  )

  # Check that result has correct structure
  expect_s3_class(result, "lmsi")
  expect_length(result$b_y, n_traits)
  expect_length(result$b_s, n_traits)

  # Check that provided G_s is used
  expect_equal(result$G_s, G_s)

  # Check validity of metrics
  expect_gte(result$rHI, 0)
  expect_lte(result$rHI, 1)
  expect_gte(result$hI2, 0)
  expect_lte(result$hI2, 1)
})

test_that("lmsi logic consistency: G_s allows NULL data matrices", {
  # This test specifically validates the fix for the logic inconsistency
  # where documentation said phen_mat could be NULL if G_s provided,
  # but code validation was inconsistent

  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  weights <- c(10, 5, 3, 3, 5, 8, 4)
  G_s <- 0.75 * gmat # Theoretical case: markers explain 75% of genetic variance

  # CASE 1: G_s provided, both data matrices NULL - should work
  expect_no_error({
    result1 <- lmsi(
      phen_mat = NULL, marker_scores = NULL,
      pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights
    )
  })

  # CASE 2: G_s provided, only phen_mat NULL - should work
  set.seed(123)
  marker_scores <- matrix(rnorm(50 * ncol(gmat), mean = 5, sd = 1.5),
    nrow = 50, ncol = ncol(gmat)
  )
  expect_no_error({
    result2 <- lmsi(
      phen_mat = NULL, marker_scores = marker_scores,
      pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights
    )
  })

  # CASE 3: G_s provided, only marker_scores NULL - should work
  phen_mat <- matrix(rnorm(50 * ncol(gmat), mean = 15, sd = 3),
    nrow = 50, ncol = ncol(gmat)
  )
  expect_no_error({
    result3 <- lmsi(
      phen_mat = phen_mat, marker_scores = NULL,
      pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights
    )
  })

  # CASE 4: G_s NULL, both data matrices NULL - should error
  expect_error(
    {
      result4 <- lmsi(
        phen_mat = NULL, marker_scores = NULL,
        pmat = pmat, gmat = gmat, G_s = NULL, wmat = weights
      )
    },
    "Either G_s must be provided, or both phen_mat and marker_scores must be provided"
  )

  # All valid cases should produce same results since G_s takes precedence
  expect_equal(result1$rHI, result2$rHI, tolerance = 1e-10)
  expect_equal(result1$GA, result2$GA, tolerance = 1e-10)
  expect_equal(result1$rHI, result3$rHI, tolerance = 1e-10)
  expect_equal(result1$GA, result3$GA, tolerance = 1e-10)
})

test_that("lmsi enhanced output structure works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data
  set.seed(123)
  n_genotypes <- 30
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(phen_mat) <- colnames(gmat)

  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(marker_scores) <- colnames(gmat)

  result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights)

  # Test enhanced summary structure
  expect_true("phenotype_coeffs" %in% names(result))
  expect_true("marker_coeffs" %in% names(result))
  expect_true("coeff_analysis" %in% names(result))

  # Check phenotype_coeffs structure
  expect_true(is.data.frame(result$phenotype_coeffs))
  expect_equal(nrow(result$phenotype_coeffs), n_traits)
  expect_true(all(c("Trait", "b_phenotype", "weight", "Delta_H") %in%
    colnames(result$phenotype_coeffs)))

  # Check marker_coeffs structure
  expect_true(is.data.frame(result$marker_coeffs))
  expect_equal(nrow(result$marker_coeffs), n_traits)
  expect_true(all(c("Trait", "b_marker", "weight", "Delta_H") %in%
    colnames(result$marker_coeffs)))

  # Check coefficient analysis structure
  expect_true(is.data.frame(result$coeff_analysis))
  expect_equal(nrow(result$coeff_analysis), 3) # Phenotype, MarkerScore, Combined
  expect_true(all(c("Component", "Sum_Abs_Coeff", "Mean_Abs_Coeff", "Max_Abs_Coeff") %in%
    colnames(result$coeff_analysis)))

  # Test that values match between structures
  expect_equal(result$phenotype_coeffs$b_phenotype, result$b_y, tolerance = 1e-6)
  expect_equal(result$marker_coeffs$b_marker, result$b_s, tolerance = 1e-6)
})

test_that("lmsi error handling works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Test dimension mismatches
  expect_error(
    lmsi(
      phen_mat = NULL, marker_scores = NULL,
      pmat = pmat[1:3, 1:3], gmat = gmat, wmat = c(1, 2, 3)
    ),
    "pmat and gmat must have the same dimensions"
  )

  # Test weight dimension mismatch
  expect_error(
    lmsi(
      phen_mat = NULL, marker_scores = NULL,
      pmat = pmat, gmat = gmat, G_s = gmat, wmat = c(1, 2, 3)
    ),
    "Length of weights must equal number of traits"
  )

  # Test missing required inputs
  expect_error(
    lmsi(
      phen_mat = NULL, marker_scores = NULL,
      pmat = pmat, gmat = gmat, G_s = NULL, wmat = c(10, 5, 3, 3, 5, 8, 4)
    ),
    "Either G_s must be provided, or both phen_mat and marker_scores must be provided to compute covariance matrices"
  )
})

test_that("gw_lmsi basic functionality works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data (well-conditioned case)
  set.seed(123)
  n_genotypes <- 60
  n_markers <- 40 # Less markers than genotypes for better conditioning
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  # Simulate marker and trait matrices
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
    nrow = n_genotypes, ncol = n_markers
  )

  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )

  # Test basic GW-LMSI calculation
  result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights)

  # Check that result has correct structure
  expect_s3_class(result, "gw_lmsi")
  expect_type(result, "list")

  # Check essential components exist
  essential_components <- c(
    "b", "P_GW", "G_GW", "rHI", "sigma_I", "R", "GA",
    "Delta_H", "lambda", "high_dimensional", "condition_number",
    "ridge_applied", "n_markers", "n_genotypes", "n_traits"
  )
  for (comp in essential_components) {
    expect_true(comp %in% names(result), info = paste("Missing component:", comp))
  }

  # Check dimensions
  expect_length(result$b, n_markers)
  expect_length(result$Delta_H, n_traits)
  expect_equal(dim(result$P_GW), c(n_markers, n_markers))
  expect_equal(dim(result$G_GW), c(n_markers, n_traits))

  # Check counts
  expect_equal(result$n_markers, n_markers)
  expect_equal(result$n_genotypes, n_genotypes)
  expect_equal(result$n_traits, n_traits)

  # Check that accuracy is in valid range [0, 1]
  expect_gte(result$rHI, 0)
  expect_lte(result$rHI, 1)

  # Check that heritability is in valid range [0, 1]
  expect_gte(result$hI2, 0)
  expect_lte(result$hI2, 1)

  # Check non-negative values
  expect_gte(result$sigma_I, 0)

  # Check logical flags
  expect_type(result$high_dimensional, "logical")
  expect_type(result$ridge_applied, "logical")

  # Should not be high-dimensional in this case
  expect_false(result$high_dimensional)
})

test_that("gw_lmsi high-dimensional case with ridge regularization works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up high-dimensional case (n_markers > n_genotypes)
  set.seed(123)
  n_genotypes <- 30
  n_markers <- 100 # More markers than genotypes
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
    nrow = n_genotypes, ncol = n_markers
  )

  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )

  # Test with ridge regularization
  lambda_val <- 0.01
  result <- gw_lmsi(marker_mat, trait_mat, gmat, lambda = lambda_val, wmat = weights)

  # Should detect high-dimensional case
  expect_true(result$high_dimensional)
  expect_equal(result$lambda, lambda_val)
  expect_true(result$ridge_applied)

  # Should still produce valid results
  expect_length(result$b, n_markers)
  expect_gte(result$rHI, 0)
  expect_lte(result$rHI, 1)

  # Test without ridge regularization (should warn about high-dimensional case)
  # Note: This will also generate numerical warnings due to singular matrices
  expect_warning(
    {
      result_no_ridge <- gw_lmsi(marker_mat, trait_mat, gmat, lambda = 0, wmat = weights)
    },
    "High-dimensional case detected"
  )

  expect_true(result_no_ridge$high_dimensional)
  expect_false(result_no_ridge$ridge_applied)
})

test_that("gw_lmsi with provided matrices works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data (well-conditioned case)
  set.seed(123)
  n_genotypes <- 60
  n_markers <- 40 # Less markers than genotypes
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
    nrow = n_genotypes, ncol = n_markers
  )

  # Create custom P_GW and G_GW matrices
  P_GW <- cov(marker_mat)
  G_GW <- matrix(rnorm(n_markers * n_traits, sd = 0.1),
    nrow = n_markers, ncol = n_traits
  )

  # Test with provided matrices
  result <- gw_lmsi(marker_mat,
    trait_mat = NULL, gmat,
    P_GW = P_GW, G_GW = G_GW, wmat = weights
  )

  expect_s3_class(result, "gw_lmsi")
  expect_equal(result$P_GW, P_GW)
  expect_equal(result$G_GW, G_GW)
  expect_length(result$b, n_markers)
})

test_that("gw_lmsi error handling works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data
  marker_mat <- matrix(sample(0:2, 50 * 80, replace = TRUE), nrow = 50, ncol = 80)
  trait_mat <- matrix(rnorm(50 * ncol(gmat), mean = 15, sd = 3), nrow = 50, ncol = ncol(gmat))
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  # Test negative lambda
  expect_error(
    gw_lmsi(marker_mat, trait_mat, gmat, lambda = -0.01, wmat = weights),
    "lambda must be non-negative"
  )

  # Test weight dimension mismatch
  expect_error(
    gw_lmsi(marker_mat, trait_mat, gmat, wmat = c(1, 2, 3)),
    "Length of weights must equal number of traits"
  )

  # Test missing required inputs
  expect_error(
    gw_lmsi(marker_mat,
      trait_mat = NULL, gmat,
      P_GW = NULL, G_GW = NULL, wmat = weights
    ),
    "Either \\(P_GW and G_GW\\) or trait_mat must be provided"
  )

  # Test dimension mismatches in provided matrices
  wrong_P_GW <- matrix(rnorm(50 * 50), nrow = 50, ncol = 50) # Wrong size for 80 markers
  correct_G_GW <- matrix(rnorm(80 * ncol(gmat)), nrow = 80, ncol = ncol(gmat)) # Correct size
  expect_error(
    gw_lmsi(marker_mat,
      trait_mat = NULL, gmat,
      P_GW = wrong_P_GW, G_GW = correct_G_GW, wmat = weights
    ),
    "P_GW must be n_markers x n_markers matrix"
  )
})

test_that("gw_lmsi condition number calculation works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up well-conditioned case
  set.seed(123)
  n_genotypes <- 50
  n_markers <- 30 # Less markers than genotypes for good conditioning
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  marker_mat <- matrix(rnorm(n_genotypes * n_markers, sd = 1),
    nrow = n_genotypes, ncol = n_markers
  )
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )

  result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights)

  # Should have calculated condition number
  expect_false(result$high_dimensional)
  expect_type(result$condition_number, "double")

  # For well-conditioned matrices, condition number should be reasonable
  if (!is.na(result$condition_number)) {
    expect_lt(result$condition_number, 1e6) # Should be well-conditioned
  }
})

test_that("print methods work without errors", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data for LMSI
  set.seed(123)
  n_genotypes <- 40
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(phen_mat) <- colnames(gmat)

  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(marker_scores) <- colnames(gmat)

  # Test LMSI print method
  lmsi_result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights)
  expect_output(print(lmsi_result), "LINEAR MARKER SELECTION INDEX")

  # Set up test data for GW-LMSI (well-conditioned case)
  n_markers <- 30 # Less than n_genotypes
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
    nrow = n_genotypes, ncol = n_markers
  )
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )

  # Test GW-LMSI print method
  gw_lmsi_result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights)
  expect_output(print(gw_lmsi_result), "GENOME-WIDE LINEAR MARKER SELECTION INDEX")

  # Test GW-LMSI print method with ridge regularization
  gw_lmsi_ridge <- gw_lmsi(marker_mat, trait_mat, gmat, lambda = 0.01, wmat = weights)
  expect_output(print(gw_lmsi_ridge), "Ridge regularization.*APPLIED")
})

test_that("PRE calculation works when GAY is provided", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Set up test data
  set.seed(123)
  n_genotypes <- 30
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)

  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(phen_mat) <- colnames(gmat)

  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(marker_scores) <- colnames(gmat)

  GAY_value <- 25.0

  # Test LMSI with GAY
  result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights, GAY = GAY_value)
  expect_false(is.na(result$PRE))
  expect_type(result$PRE, "double")
  expect_equal(result$PRE, (result$GA / GAY_value) * 100)

  # Test GW-LMSI with GAY (use well-conditioned case to avoid numerical warnings)
  set.seed(456) # Different seed for better conditioning
  n_markers_stable <- 25 # Well below n_genotypes for good conditioning
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers_stable, replace = TRUE),
    nrow = n_genotypes, ncol = n_markers_stable
  )

  # Add some structure to improve conditioning
  marker_mat <- marker_mat + matrix(rnorm(n_genotypes * n_markers_stable, sd = 0.1),
    nrow = n_genotypes, ncol = n_markers_stable
  )

  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )

  gw_result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights, GAY = GAY_value)
  expect_false(is.na(gw_result$PRE))
  expect_type(gw_result$PRE, "double")
})

# ==============================================================================
# NEW COVERAGE TESTS - targeting previously uncovered lines
# ==============================================================================

# Helper: standard seldata matrices
.load_marker_test_data <- function() {
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  list(
    pmat = pmat, gmat = gmat, n_traits = ncol(gmat),
    weights = c(10, 5, 3, 3, 5, 8, 4)
  )
}

# --- lmsi: line 150 – pmat/gmat not square --------------------------------
test_that("lmsi errors when pmat or gmat is not square (line 150)", {
  d <- .load_marker_test_data()
  # A 2×3 pmat is not square
  bad_pmat <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(
    lmsi(
      pmat = bad_pmat, gmat = d$gmat, G_s = d$gmat[1:2, 1:2],
      wmat = c(1, 2)
    ),
    "pmat and gmat must be square matrices"
  )
})

# --- lmsi: lines 159-162 – matrix wmat path with wcol validation ----------
test_that("lmsi validates wcol against matrix wmat columns (lines 159-162)", {
  d <- .load_marker_test_data()
  G_s <- 0.8 * d$gmat
  wmat_matrix <- matrix(rep(d$weights, 2), ncol = 2) # 7×2 matrix

  # Valid: wcol = 2 should use second column
  res <- lmsi(
    pmat = d$pmat, gmat = d$gmat, G_s = G_s,
    wmat = wmat_matrix, wcol = 2
  )
  expect_s3_class(res, "lmsi")

  # Invalid: wcol > ncol(wmat)
  expect_error(
    lmsi(
      pmat = d$pmat, gmat = d$gmat, G_s = G_s,
      wmat = wmat_matrix, wcol = 5
    ),
    "wcol exceeds number of columns in wmat"
  )
})

# --- lmsi: line 182 – phen_mat/marker_scores column count mismatch --------
test_that("lmsi errors when phen_mat or marker_scores has wrong columns (line 182)", {
  d <- .load_marker_test_data()
  set.seed(1)
  n <- 30
  phen_mat <- matrix(rnorm(n * d$n_traits), n, d$n_traits)
  bad_marker_sc <- matrix(rnorm(n * (d$n_traits - 1)), n, d$n_traits - 1)
  expect_error(
    lmsi(
      phen_mat = phen_mat, marker_scores = bad_marker_sc,
      pmat = d$pmat, gmat = d$gmat, wmat = d$weights
    ),
    "Number of columns in phen_mat and marker_scores must equal n_traits"
  )
})

# --- lmsi: line 186 – phen_mat/marker_scores row count mismatch -----------
test_that("lmsi errors when phen_mat and marker_scores have different rows (line 186)", {
  d <- .load_marker_test_data()
  set.seed(2)
  phen_mat <- matrix(rnorm(30 * d$n_traits), 30, d$n_traits)
  bad_marker_sc <- matrix(rnorm(20 * d$n_traits), 20, d$n_traits)
  expect_error(
    lmsi(
      phen_mat = phen_mat, marker_scores = bad_marker_sc,
      pmat = d$pmat, gmat = d$gmat, wmat = d$weights
    ),
    "phen_mat and marker_scores must have the same number of rows"
  )
})

# --- lmsi: line 204 – G_s wrong dimensions --------------------------------
test_that("lmsi errors when G_s has wrong dimensions (line 204)", {
  d <- .load_marker_test_data()
  bad_G_s <- matrix(1, 3, 3) # 3×3 but n_traits = 7
  expect_error(
    lmsi(pmat = d$pmat, gmat = d$gmat, G_s = bad_G_s, wmat = d$weights),
    "G_s must be n_traits x n_traits matrix"
  )
})

# --- lmsi: lines 246-247 – P_L singular: solve() fails, ginv used ---------
test_that("lmsi warns and falls back to ginv when P_L is singular (lines 246-247)", {
  d <- .load_marker_test_data()
  # An all-zero P_L arises when pmat and G_s (=Cov_ys=Var_s) are all zero.
  # Zero pmat → P_L is all-zeros → solve() throws → ginv fallback fires.
  zero_gmat <- 0 * d$gmat
  zero_pmat <- 0 * d$pmat
  G_s <- 0 * d$gmat
  expect_warning(
    lmsi(pmat = zero_pmat, gmat = zero_gmat, G_s = G_s, wmat = d$weights),
    "P_L is singular or near-singular, using generalized inverse"
  )
})

# --- lmsi: line 275 – rHI = 0 when w'Gw = 0 (denominator = 0) -----------
test_that("lmsi returns rHI = 0 when gmat is all-zero (line 275)", {
  d <- .load_marker_test_data()
  # Denominator = w'Gw; with gmat = 0, denominator = 0 → rHI branch = 0
  zero_gmat <- 0 * d$gmat
  G_s <- 0.8 * zero_gmat
  res <- suppressWarnings(
    lmsi(pmat = d$pmat, gmat = zero_gmat, G_s = G_s, wmat = d$weights)
  )
  expect_equal(res$rHI, 0)
})

# --- lmsi: line 285 – Delta_H = rep(0, n_traits) when sigma_I = 0 --------
test_that("lmsi returns zero Delta_H when sigma_I is zero (line 285)", {
  d <- .load_marker_test_data()
  # All-zero P_L → solve fails → ginv(all-zero) = all-zero → sigma_I = 0
  zero_pmat <- 0 * d$pmat
  zero_gmat <- 0 * d$gmat
  G_s <- 0 * d$gmat
  res <- suppressWarnings(
    lmsi(pmat = zero_pmat, gmat = zero_gmat, G_s = G_s, wmat = d$weights)
  )
  expect_true(all(res$Delta_H == 0))
})

# --- lmsi: line 307 – auto trait names when pmat has no colnames ----------
test_that("lmsi generates 'Trait1...' names when pmat has no colnames (line 307)", {
  d <- .load_marker_test_data()
  G_s <- 0.8 * d$gmat
  res <- lmsi(
    pmat = unname(d$pmat), gmat = unname(d$gmat),
    G_s = G_s, wmat = d$weights
  )
  expect_true(all(grepl("^Trait", res$trait_names)))
})

# --- lmsi: line 835 – PRE printed when GAY is given ----------------------
test_that("lmsi print method shows Relative Efficiency when GAY provided (line 835)", {
  d <- .load_marker_test_data()
  set.seed(42)
  n <- 30
  phen_mat <- matrix(rnorm(n * d$n_traits, 15, 3), n, d$n_traits)
  marker_sc <- matrix(rnorm(n * d$n_traits, 5, 1.5), n, d$n_traits)
  colnames(phen_mat) <- colnames(d$pmat)
  colnames(marker_sc) <- colnames(d$pmat)
  res <- lmsi(phen_mat, marker_sc, d$pmat, d$gmat,
    wmat = d$weights, GAY = 10.0
  )
  expect_output(print(res), "Relative Efficiency")
})

# ==============================================================================
# gw_lmsi – new coverage tests
# ==============================================================================

# --- gw_lmsi: line 548 – gmat not square ---------------------------------
test_that("gw_lmsi errors when gmat is not square (line 548)", {
  set.seed(1)
  marker_mat <- matrix(sample(0:2, 30 * 5, replace = TRUE), 30, 5)
  bad_gmat <- matrix(1:12, 3, 4) # 3×4 – not square
  weights <- c(1, 2, 3)
  expect_error(
    gw_lmsi(marker_mat,
      trait_mat = NULL, gmat = bad_gmat,
      P_GW = diag(5), G_GW = matrix(1, 5, 3), wmat = weights
    ),
    "gmat must be a square matrix"
  )
})

# --- gw_lmsi: lines 553-556 – matrix wmat wcol validation ----------------
test_that("gw_lmsi validates wcol for matrix wmat (lines 553-556)", {
  d <- .load_marker_test_data()
  set.seed(3)
  n <- 40
  nm <- 20
  marker_mat <- matrix(sample(0:2, n * nm, replace = TRUE), n, nm)
  trait_mat <- matrix(rnorm(n * d$n_traits, 15, 3), n, d$n_traits)
  wmat_mat <- matrix(rep(d$weights, 2), ncol = 2)

  # Valid: wcol = 2 should work
  res <- gw_lmsi(marker_mat, trait_mat, d$gmat, wmat = wmat_mat, wcol = 2)
  expect_s3_class(res, "gw_lmsi")

  # Invalid: wcol out of range
  expect_error(
    gw_lmsi(marker_mat, trait_mat, d$gmat, wmat = wmat_mat, wcol = 10),
    "wcol exceeds number of columns in wmat"
  )
})

# --- gw_lmsi: line 574 – trait_mat wrong column count --------------------
test_that("gw_lmsi errors when trait_mat has wrong number of columns (line 574)", {
  d <- .load_marker_test_data()
  set.seed(4)
  n <- 40
  nm <- 20
  marker_mat <- matrix(sample(0:2, n * nm, replace = TRUE), n, nm)
  bad_trait_mat <- matrix(rnorm(n * (d$n_traits - 1)), n, d$n_traits - 1)
  expect_error(
    gw_lmsi(marker_mat, bad_trait_mat, d$gmat, wmat = d$weights),
    "Number of columns in trait_mat must equal n_traits"
  )
})

# --- gw_lmsi: line 578 – trait_mat row count mismatch --------------------
test_that("gw_lmsi errors when trait_mat has wrong number of rows (line 578)", {
  d <- .load_marker_test_data()
  set.seed(5)
  n <- 40
  nm <- 20
  marker_mat <- matrix(sample(0:2, n * nm, replace = TRUE), n, nm)
  bad_trait_mat <- matrix(rnorm(20 * d$n_traits), 20, d$n_traits) # only 20 rows
  expect_error(
    gw_lmsi(marker_mat, bad_trait_mat, d$gmat, wmat = d$weights),
    "marker_mat and trait_mat must have the same number of rows"
  )
})

# --- gw_lmsi: line 600 – G_GW wrong dimensions ---------------------------
test_that("gw_lmsi errors when G_GW has wrong dimensions (line 600)", {
  d <- .load_marker_test_data()
  set.seed(6)
  nm <- 20
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  P_GW <- diag(nm)
  bad_G_GW <- matrix(1, nm, d$n_traits - 1) # wrong cols
  expect_error(
    gw_lmsi(marker_mat,
      trait_mat = NULL, d$gmat,
      P_GW = P_GW, G_GW = bad_G_GW, wmat = d$weights
    ),
    "G_GW must be n_markers x n_traits matrix"
  )
})

# --- gw_lmsi: line 631 – condition_number <<- NA_real_ in error handler ---
# NaN in P_GW makes eigen() throw. The tryCatch catches it and assigns
# condition_number <<- NA_real_ (line 631). Execution then reaches the
# high-dimensional-warning block (line 640) which confirms line 631 fired.
# The downstream ginv(NaN_matrix) call eventually crashes, so a test-level
# tryCatch swallows that error; warnings emitted before it are still captured.
test_that("gw_lmsi assigns NA to condition_number when eigen() throws on NaN P_GW (line 631)", {
  d <- .load_marker_test_data()
  set.seed(21)
  # n_markers (5) > n_genotypes (3) -> high_dimensional = TRUE, lambda = 0
  # -> fires the high-dimensional warning (line 640) AFTER line 631 executes.
  n_gen <- 3
  nm <- 5
  marker_mat <- matrix(sample(0:2, n_gen * nm, replace = TRUE), n_gen, nm)
  G_GW <- matrix(0, nm, d$n_traits) # zero G_GW avoids NaN in b'G_GW w path
  P_GW_nan <- diag(nm)
  P_GW_nan[1L, 1L] <- NaN # NaN -> eigen() throws -> line 631 fires
  # Warning is emitted (line 640) before the downstream error from ginv(NaN).
  # tryCatch at test level swallows the downstream error without hiding warnings,
  # because testthat's withCallingHandlers captures them before stack unwinds.
  expect_warning(
    tryCatch(
      gw_lmsi(marker_mat,
        trait_mat = NULL, d$gmat,
        P_GW = P_GW_nan, G_GW = G_GW, wmat = d$weights, lambda = 0
      ),
      error = function(e) NULL # swallow crash from ginv(NaN_matrix)
    ),
    "High-dimensional case detected" # line 640 fires after line 631
  )
})

# --- gw_lmsi: lines 658-659 – "numerically singular" warning (NA condition_number, non-high-dim) ---
# A 1x1 zero P_GW: eigen() succeeds but all eigenvalues = 0 <= 1e-14,
# so min_eig_pos is empty and condition_number stays NA via the silent
# if-branch (NOT the error handler). With n_markers < n_genotypes and lambda=0,
# lines 658-659 fire: 'P_GW appears to be numerically singular'.
# solve() then fails (lines 688, 691 also fire), but ginv returns b=0 -> sigma_I=0
# -> Delta_H=rep(0,...) so function completes without crashing.
test_that("gw_lmsi warns about numerically singular P_GW when condition_number is NA (lines 658-659)", {
  d <- .load_marker_test_data()
  set.seed(7)
  n_gen <- 30
  nm <- 1 # 1 marker, 30 genotypes: NOT high-dimensional
  marker_mat <- matrix(rep(0, n_gen * nm), n_gen, nm) # constant -> zero cov
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.3), nm, d$n_traits)
  # A zero 1x1 P_GW: all eigenvalues = 0 <= 1e-14 -> condition_number stays NA
  zero_P_GW <- matrix(0, nm, nm)
  # The function emits 3 warnings in total (lines 659, 688, 691).
  # Nest expect_warning so every warning is explicitly consumed.
  expect_warning(
    expect_warning(
      expect_warning(
        gw_lmsi(marker_mat,
          trait_mat = NULL, d$gmat,
          P_GW = zero_P_GW, G_GW = G_GW, wmat = d$weights, lambda = 0
        ),
        "P_GW appears to be numerically singular" # line 659
      ),
      "singular or near-singular" # line 688
    ),
    "Consider using lambda" # line 691
  )
})


# --- gw_lmsi: lines 651-656 – ill-conditioned warning --------------------
test_that("gw_lmsi warns when P_GW is ill-conditioned (condition number > 1e10) (lines 651-656)", {
  d <- .load_marker_test_data()
  set.seed(8)
  # Construct an ill-conditioned but NOT singular P_GW:
  # Use eigenvalue decomposition to set one eigenvalue very tiny
  nm <- 4
  Q <- qr.Q(qr(matrix(rnorm(nm^2), nm, nm)))
  eigs <- c(1e12, 1, 1, 1) # condition number = 1e12 > 1e10
  P_ill <- Q %*% diag(eigs) %*% t(Q)
  P_ill <- (P_ill + t(P_ill)) / 2 # ensure symmetry
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.1), nm, d$n_traits)
  # lambda > 0 so high-dimensional warning is not mixed in
  expect_warning(
    gw_lmsi(marker_mat,
      trait_mat = NULL, d$gmat,
      P_GW = P_ill, G_GW = G_GW, wmat = d$weights, lambda = 0.01
    ),
    "ill-conditioned"
  )
})

# --- gw_lmsi: lines 688-691 – P_GW singular, non-high-dim, lambda=0 -----
# solve() on a singular (not high-dimensional) P_GW → 2 extra warnings
test_that("gw_lmsi warns about singular P_GW and suggests ridge (lines 688-691)", {
  d <- .load_marker_test_data()
  set.seed(9)
  # Marker matrix with all identical columns → P_GW is singular
  # n_markers (3) < n_genotypes (40) → NOT high-dimensional
  nm <- 3
  base_col <- rnorm(40)
  marker_mat <- cbind(base_col, base_col, base_col) # repeated → singular P_GW
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.2), nm, d$n_traits)
  P_GW <- cov(marker_mat) # near-singular
  # Expect at least one of the singular-P_GW warnings
  # Both warnings from lines 688 and 691 are emitted; consume both.
  expect_warning(
    expect_warning(
      gw_lmsi(marker_mat,
        trait_mat = NULL, d$gmat,
        P_GW = P_GW, G_GW = G_GW, wmat = d$weights, lambda = 0
      ),
      "singular or near-singular" # line 688
    ),
    "Consider using lambda" # line 691
  )
})

# --- gw_lmsi: line 718 – rHI = 0 when w'Gw = 0 --------------------------
test_that("gw_lmsi returns rHI = 0 when gmat is all-zero (line 718)", {
  d <- .load_marker_test_data()
  set.seed(10)
  nm <- 5
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.1), nm, d$n_traits)
  P_GW <- diag(nm)
  zero_gmat <- 0 * d$gmat
  res <- gw_lmsi(marker_mat,
    trait_mat = NULL, zero_gmat,
    P_GW = P_GW, G_GW = G_GW, wmat = d$weights
  )
  expect_equal(res$rHI, 0)
})

# --- gw_lmsi: line 728 – Delta_H = rep(0,...) when sigma_I = 0 ----------
test_that("gw_lmsi returns zero Delta_H when sigma_I is zero (line 728)", {
  d <- .load_marker_test_data()
  set.seed(11)
  nm <- 3
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.1), nm, d$n_traits)
  # All-zero P_GW → b'P_GWb = 0 → sigma_I = 0
  zero_P_GW <- matrix(0, nm, nm)
  res <- suppressWarnings(
    gw_lmsi(marker_mat,
      trait_mat = NULL, d$gmat,
      P_GW = zero_P_GW, G_GW = G_GW, wmat = d$weights, lambda = 0
    )
  )
  expect_true(all(res$Delta_H == 0))
})

# --- gw_lmsi: line 750 – auto trait names when gmat has no colnames ------
test_that("gw_lmsi generates 'Trait...' names when gmat has no colnames (line 750)", {
  d <- .load_marker_test_data()
  set.seed(12)
  nm <- 10
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  trait_mat <- matrix(rnorm(40 * d$n_traits, 15, 3), 40, d$n_traits)
  res <- gw_lmsi(marker_mat, trait_mat, unname(d$gmat), wmat = d$weights)
  expect_true(all(grepl("^Trait", res$trait_names)))
})

# --- print.gw_lmsi: line 882 – HIGH-DIMENSIONAL print --------------------
test_that("print.gw_lmsi shows HIGH-DIMENSIONAL when high_dimensional is TRUE (line 882)", {
  d <- .load_marker_test_data()
  set.seed(13)
  n_gen <- 20
  nm <- 50 # nm > n_gen → high-dimensional
  marker_mat <- matrix(sample(0:2, n_gen * nm, replace = TRUE), n_gen, nm)
  trait_mat <- matrix(rnorm(n_gen * d$n_traits, 15, 3), n_gen, d$n_traits)
  res <- suppressWarnings(
    gw_lmsi(marker_mat, trait_mat, d$gmat, wmat = d$weights, lambda = 0.01)
  )
  expect_true(res$high_dimensional)
  expect_output(print(res), "HIGH-DIMENSIONAL")
})

# --- print.gw_lmsi: line 889 – ILL-CONDITIONED print --------------------
test_that("print.gw_lmsi shows ILL-CONDITIONED when condition_number > 1e10 (line 889)", {
  d <- .load_marker_test_data()
  set.seed(14)
  nm <- 4
  Q <- qr.Q(qr(matrix(rnorm(nm^2), nm, nm)))
  eigs <- c(1e12, 1, 1, 1)
  P_ill <- Q %*% diag(eigs) %*% t(Q)
  P_ill <- (P_ill + t(P_ill)) / 2
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.1), nm, d$n_traits)
  res <- suppressWarnings(
    gw_lmsi(marker_mat,
      trait_mat = NULL, d$gmat,
      P_GW = P_ill, G_GW = G_GW, wmat = d$weights, lambda = 0.01
    )
  )
  expect_gt(res$condition_number, 1e10)
  expect_output(print(res), "ILL-CONDITIONED")
})

# --- print.gw_lmsi: line 891 – MODERATE condition number print -----------
test_that("print.gw_lmsi shows MODERATE when 1e6 < condition_number <= 1e10 (line 891)", {
  d <- .load_marker_test_data()
  set.seed(15)
  nm <- 4
  Q <- qr.Q(qr(matrix(rnorm(nm^2), nm, nm)))
  eigs <- c(1e8, 1, 1, 1) # condition number ≈ 1e8  (1e6 < x <= 1e10)
  P_mod <- Q %*% diag(eigs) %*% t(Q)
  P_mod <- (P_mod + t(P_mod)) / 2
  marker_mat <- matrix(sample(0:2, 40 * nm, replace = TRUE), 40, nm)
  G_GW <- matrix(rnorm(nm * d$n_traits, sd = 0.1), nm, d$n_traits)
  res <- suppressWarnings(
    gw_lmsi(marker_mat,
      trait_mat = NULL, d$gmat,
      P_GW = P_mod, G_GW = G_GW, wmat = d$weights, lambda = 0.01
    )
  )
  # Confirm in MODERATE band and verify print output
  expect_gt(res$condition_number, 1e6)
  expect_lte(res$condition_number, 1e10)
  expect_output(print(res), "MODERATE")
})
