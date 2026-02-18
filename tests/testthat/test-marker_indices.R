test_that("lmsi basic functionality works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data
  set.seed(123)
  n_genotypes <- 50
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  # Simulate phenotype and marker score matrices
  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                     nrow = n_genotypes, ncol = n_traits)
  colnames(phen_mat) <- colnames(gmat)
  
  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
                          nrow = n_genotypes, ncol = n_traits)
  colnames(marker_scores) <- colnames(gmat)
  
  # Test basic LMSI calculation
  result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights)
  
  # Check that result has correct structure
  expect_s3_class(result, "lmsi")
  expect_type(result, "list")
  
  # Check essential components exist
  essential_components <- c("b_y", "b_s", "b_combined", "P_L", "G_L", "G_s", 
                           "rHI", "sigma_I", "R", "GA", "Delta_H", "summary",
                           "phenotype_coeffs", "marker_coeffs", "coeff_analysis")
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
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  # Create custom G_s matrix
  G_s <- 0.8 * gmat  # Assume markers explain 80% of genetic variance
  
  # Test LMSI with provided G_s (no data matrices needed)
  result <- lmsi(phen_mat = NULL, marker_scores = NULL, 
                 pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights)
  
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
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  G_s <- 0.75 * gmat  # Theoretical case: markers explain 75% of genetic variance
  
  # CASE 1: G_s provided, both data matrices NULL - should work
  expect_no_error({
    result1 <- lmsi(phen_mat = NULL, marker_scores = NULL,
                    pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights)
  })
  
  # CASE 2: G_s provided, only phen_mat NULL - should work  
  set.seed(123)
  marker_scores <- matrix(rnorm(50 * ncol(gmat), mean = 5, sd = 1.5),
                          nrow = 50, ncol = ncol(gmat))
  expect_no_error({
    result2 <- lmsi(phen_mat = NULL, marker_scores = marker_scores,
                    pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights)
  })
  
  # CASE 3: G_s provided, only marker_scores NULL - should work
  phen_mat <- matrix(rnorm(50 * ncol(gmat), mean = 15, sd = 3),
                     nrow = 50, ncol = ncol(gmat))
  expect_no_error({
    result3 <- lmsi(phen_mat = phen_mat, marker_scores = NULL,
                    pmat = pmat, gmat = gmat, G_s = G_s, wmat = weights)
  })
  
  # CASE 4: G_s NULL, both data matrices NULL - should error
  expect_error({
    result4 <- lmsi(phen_mat = NULL, marker_scores = NULL,
                    pmat = pmat, gmat = gmat, G_s = NULL, wmat = weights)
  }, "Either G_s must be provided, or both phen_mat and marker_scores must be provided")
  
  # All valid cases should produce same results since G_s takes precedence
  expect_equal(result1$rHI, result2$rHI, tolerance = 1e-10)
  expect_equal(result1$GA, result2$GA, tolerance = 1e-10)
  expect_equal(result1$rHI, result3$rHI, tolerance = 1e-10)
  expect_equal(result1$GA, result3$GA, tolerance = 1e-10)
})

test_that("lmsi enhanced output structure works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data
  set.seed(123)
  n_genotypes <- 30
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                     nrow = n_genotypes, ncol = n_traits)
  colnames(phen_mat) <- colnames(gmat)
  
  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
                          nrow = n_genotypes, ncol = n_traits)
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
  expect_equal(nrow(result$coeff_analysis), 3)  # Phenotype, MarkerScore, Combined
  expect_true(all(c("Component", "Sum_Abs_Coeff", "Mean_Abs_Coeff", "Max_Abs_Coeff") %in% 
                  colnames(result$coeff_analysis)))
  
  # Test that values match between structures
  expect_equal(result$phenotype_coeffs$b_phenotype, result$b_y, tolerance = 1e-6)
  expect_equal(result$marker_coeffs$b_marker, result$b_s, tolerance = 1e-6)
})

test_that("lmsi error handling works", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Test dimension mismatches
  expect_error(lmsi(phen_mat = NULL, marker_scores = NULL, 
                    pmat = pmat[1:3, 1:3], gmat = gmat, wmat = c(1, 2, 3)),
               "pmat and gmat must have the same dimensions")
  
  # Test weight dimension mismatch
  expect_error(lmsi(phen_mat = NULL, marker_scores = NULL, 
                    pmat = pmat, gmat = gmat, G_s = gmat, wmat = c(1, 2, 3)),
               "Length of weights must equal number of traits")
  
  # Test missing required inputs
  expect_error(lmsi(phen_mat = NULL, marker_scores = NULL, 
                    pmat = pmat, gmat = gmat, G_s = NULL, wmat = c(10, 5, 3, 3, 5, 8, 4)),
               "Either G_s must be provided, or both phen_mat and marker_scores must be provided to compute covariance matrices")
})

test_that("gw_lmsi basic functionality works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data (well-conditioned case)
  set.seed(123)
  n_genotypes <- 60
  n_markers <- 40   # Less markers than genotypes for better conditioning
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  # Simulate marker and trait matrices
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
                       nrow = n_genotypes, ncol = n_markers)
  
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                      nrow = n_genotypes, ncol = n_traits)
  
  # Test basic GW-LMSI calculation
  result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights)
  
  # Check that result has correct structure
  expect_s3_class(result, "gw_lmsi")
  expect_type(result, "list")
  
  # Check essential components exist
  essential_components <- c("b", "P_GW", "G_GW", "rHI", "sigma_I", "R", "GA", 
                           "Delta_H", "lambda", "high_dimensional", "condition_number",
                           "ridge_applied", "n_markers", "n_genotypes", "n_traits")
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
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up high-dimensional case (n_markers > n_genotypes)
  set.seed(123)
  n_genotypes <- 30
  n_markers <- 100    # More markers than genotypes
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
                       nrow = n_genotypes, ncol = n_markers)
  
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                      nrow = n_genotypes, ncol = n_traits)
  
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
  expect_warning({
    result_no_ridge <- gw_lmsi(marker_mat, trait_mat, gmat, lambda = 0, wmat = weights)
  }, "High-dimensional case detected")
  
  expect_true(result_no_ridge$high_dimensional)
  expect_false(result_no_ridge$ridge_applied)
})

test_that("gw_lmsi with provided matrices works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data (well-conditioned case)
  set.seed(123)
  n_genotypes <- 60
  n_markers <- 40   # Less markers than genotypes
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
                       nrow = n_genotypes, ncol = n_markers)
  
  # Create custom P_GW and G_GW matrices
  P_GW <- cov(marker_mat)
  G_GW <- matrix(rnorm(n_markers * n_traits, sd = 0.1), 
                 nrow = n_markers, ncol = n_traits)
  
  # Test with provided matrices
  result <- gw_lmsi(marker_mat, trait_mat = NULL, gmat, 
                    P_GW = P_GW, G_GW = G_GW, wmat = weights)
  
  expect_s3_class(result, "gw_lmsi")
  expect_equal(result$P_GW, P_GW)
  expect_equal(result$G_GW, G_GW)
  expect_length(result$b, n_markers)
})

test_that("gw_lmsi error handling works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data
  marker_mat <- matrix(sample(0:2, 50 * 80, replace = TRUE), nrow = 50, ncol = 80)
  trait_mat <- matrix(rnorm(50 * ncol(gmat), mean = 15, sd = 3), nrow = 50, ncol = ncol(gmat))
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  # Test negative lambda
  expect_error(gw_lmsi(marker_mat, trait_mat, gmat, lambda = -0.01, wmat = weights),
               "lambda must be non-negative")
  
  # Test weight dimension mismatch
  expect_error(gw_lmsi(marker_mat, trait_mat, gmat, wmat = c(1, 2, 3)),
               "Length of weights must equal number of traits")
  
  # Test missing required inputs
  expect_error(gw_lmsi(marker_mat, trait_mat = NULL, gmat, 
                       P_GW = NULL, G_GW = NULL, wmat = weights),
               "Either \\(P_GW and G_GW\\) or trait_mat must be provided")
  
  # Test dimension mismatches in provided matrices
  wrong_P_GW <- matrix(rnorm(50 * 50), nrow = 50, ncol = 50)  # Wrong size for 80 markers
  correct_G_GW <- matrix(rnorm(80 * ncol(gmat)), nrow = 80, ncol = ncol(gmat))  # Correct size
  expect_error(gw_lmsi(marker_mat, trait_mat = NULL, gmat, 
                       P_GW = wrong_P_GW, G_GW = correct_G_GW, wmat = weights),
               "P_GW must be n_markers x n_markers matrix")
})

test_that("gw_lmsi condition number calculation works", {
  # Load test data
  data("seldata", package = "selection.index")
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up well-conditioned case
  set.seed(123)
  n_genotypes <- 50
  n_markers <- 30    # Less markers than genotypes for good conditioning
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  marker_mat <- matrix(rnorm(n_genotypes * n_markers, sd = 1), 
                       nrow = n_genotypes, ncol = n_markers)
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                      nrow = n_genotypes, ncol = n_traits)
  
  result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights)
  
  # Should have calculated condition number
  expect_false(result$high_dimensional)
  expect_type(result$condition_number, "double")
  
  # For well-conditioned matrices, condition number should be reasonable
  if (!is.na(result$condition_number)) {
    expect_lt(result$condition_number, 1e6)  # Should be well-conditioned
  }
})

test_that("print methods work without errors", {
  # Load test data
  data("seldata", package = "selection.index")
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data for LMSI
  set.seed(123)
  n_genotypes <- 40
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                     nrow = n_genotypes, ncol = n_traits)
  colnames(phen_mat) <- colnames(gmat)
  
  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
                          nrow = n_genotypes, ncol = n_traits)
  colnames(marker_scores) <- colnames(gmat)
  
  # Test LMSI print method
  lmsi_result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights)
  expect_output(print(lmsi_result), "LINEAR MARKER SELECTION INDEX")
  
  # Set up test data for GW-LMSI (well-conditioned case)
  n_markers <- 30    # Less than n_genotypes
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
                       nrow = n_genotypes, ncol = n_markers)
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                      nrow = n_genotypes, ncol = n_traits)
  
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
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Set up test data
  set.seed(123)
  n_genotypes <- 30
  n_traits <- ncol(gmat)
  weights <- c(10, 5, 3, 3, 5, 8, 4)
  
  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                     nrow = n_genotypes, ncol = n_traits)
  colnames(phen_mat) <- colnames(gmat)
  
  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
                          nrow = n_genotypes, ncol = n_traits)
  colnames(marker_scores) <- colnames(gmat)
  
  GAY_value <- 25.0
  
  # Test LMSI with GAY
  result <- lmsi(phen_mat, marker_scores, pmat, gmat, wmat = weights, GAY = GAY_value)
  expect_false(is.na(result$PRE))
  expect_type(result$PRE, "double")
  expect_equal(result$PRE, (result$GA / GAY_value) * 100)
  
  # Test GW-LMSI with GAY (use well-conditioned case to avoid numerical warnings)
  set.seed(456)  # Different seed for better conditioning
  n_markers_stable <- 25  # Well below n_genotypes for good conditioning
  marker_mat <- matrix(sample(0:2, n_genotypes * n_markers_stable, replace = TRUE),
                       nrow = n_genotypes, ncol = n_markers_stable)
  
  # Add some structure to improve conditioning
  marker_mat <- marker_mat + matrix(rnorm(n_genotypes * n_markers_stable, sd = 0.1),
                                    nrow = n_genotypes, ncol = n_markers_stable)
  
  trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                      nrow = n_genotypes, ncol = n_traits)
  
  gw_result <- gw_lmsi(marker_mat, trait_mat, gmat, wmat = weights, GAY = GAY_value)
  expect_false(is.na(gw_result$PRE))
  expect_type(gw_result$PRE, "double")
})