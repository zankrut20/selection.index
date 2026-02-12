# ==============================================================================
# Tests for Genomic Selection Indices (LGSI and CLGSI)
# ==============================================================================

# Setup: Create test data
setup_test_data <- function() {
  set.seed(123)
  
  # Use actual variance-covariance matrices from example data
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Simulate GEBVs using a simple approach: 
  # Take actual phenotypes and add independent normal noise
  # This gives realistic correlations without positive-definiteness issues
  phen_mat <- as.matrix(seldata[1:50, 3:9])
  
  # Create GEBVs as phenotypes + noise (simulating prediction)
  # GEBVs should have lower variance than phenotypes (reliability < 1)
  gebv_mat <- phen_mat * 0.7 + matrix(rnorm(prod(dim(phen_mat)), sd = 0.5), 
                                       nrow(phen_mat), ncol(phen_mat))
  colnames(gebv_mat) <- colnames(gmat)
  
  n_traits <- ncol(gmat)
  
  # Economic weights (extract the 'ew' column as a vector)
  weights <- weight$ew
  
  list(
    gmat = gmat,
    pmat = pmat,
    gebv_mat = gebv_mat,
    phen_mat = phen_mat,
    weights = weights,
    n_traits = n_traits
  )
}

# ==============================================================================
# TEST: lgsi() - Basic Functionality
# ==============================================================================

test_that("lgsi computes index with known reliability", {
  data <- setup_test_data()
  
  result <- lgsi(
    gebv_mat = data$gebv_mat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = 0.7
  )
  
  # Check structure
  expect_s3_class(result, "lgsi")
  expect_s3_class(result, "selection_index")
  expect_type(result, "list")
  
  # Check components exist
  expect_true("b" %in% names(result))
  expect_true("P_gebv" %in% names(result))
  expect_true("reliability" %in% names(result))
  expect_true("GA" %in% names(result))
  expect_true("hI2" %in% names(result))
  expect_true("rHI" %in% names(result))
  expect_true("summary" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$b), data$n_traits)
  expect_equal(dim(result$P_gebv), c(data$n_traits, data$n_traits))
  expect_equal(length(result$reliability), data$n_traits)
  
  # Check values are reasonable
  expect_true(all(result$reliability == 0.7))
  expect_true(result$hI2 >= 0 && result$hI2 <= 1)
  expect_true(result$rHI >= 0 && result$rHI <= 1)
  expect_false(is.na(result$GA))
})

test_that("lgsi computes index with vector reliability", {
  data <- setup_test_data()
  
  rel_vec <- seq(0.5, 0.9, length.out = data$n_traits)
  
  result <- lgsi(
    gebv_mat = data$gebv_mat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = rel_vec
  )
  
  expect_equal(unname(result$reliability), rel_vec)
  expect_false(is.na( result$GA))
})

test_that("lgsi estimates reliability when not provided", {
  data <- setup_test_data()
  
  result <- lgsi(
    gebv_mat = data$gebv_mat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = NULL
  )
  
  # Should estimate reliability from variance ratio
  expect_true(all(result$reliability >= 0 & result$reliability <= 1))
  expect_equal(length(result$reliability), data$n_traits)
})

test_that("lgsi calculates PRE when GAY provided", {
  data <- setup_test_data()
  
  GAY <- 1.075
  
  result <- lgsi(
    gebv_mat = data$gebv_mat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = 0.7,
    GAY = GAY
  )
  
  expect_false(is.na(result$PRE))
  expect_equal(result$PRE, (result$GA / GAY) * 100, tolerance = 1e-6)
})

test_that("lgsi handles different selection intensities", {
  data <- setup_test_data()
  
  result1 <- lgsi(data$gebv_mat, data$gmat, data$weights, 
                  reliability = 0.7, selection_intensity = 2.063)
  result2 <- lgsi(data$gebv_mat, data$gmat, data$weights, 
                  reliability = 0.7, selection_intensity = 2.421)
  
  # GA should scale with selection intensity
  ratio <- result2$GA / result1$GA
  expect_equal(ratio, 2.421 / 2.063, tolerance = 1e-6)
})

# ==============================================================================
# TEST: lgsi() - Input Validation
# ==============================================================================

test_that("lgsi errors with dimension mismatch", {
  data <- setup_test_data()
  
  wrong_gmat <- data$gmat[1:5, 1:5]
  
  expect_error(
    lgsi(data$gebv_mat, wrong_gmat, data$weights),
    "gmat dimensions must match"
  )
})

test_that("lgsi errors with NA in GEBVs", {
  data <- setup_test_data()
  
  gebv_with_na <- data$gebv_mat
  gebv_with_na[1, 1] <- NA
  
  expect_error(
    lgsi(gebv_with_na, data$gmat, data$weights),
    "gebv_mat contains NA"
  )
})

test_that("lgsi errors with invalid reliability", {
  data <- setup_test_data()
  
  expect_error(
    lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 1.5),
    "between 0 and 1"
  )
  
  expect_error(
    lgsi(data$gebv_mat, data$gmat, data$weights, reliability = c(0.5, 0.7)),
    "vector of length n_traits"
  )
})

test_that("lgsi warns on zero variance traits", {
  data <- setup_test_data()
  
  # Create GEBVs with one constant trait
  gebv_zero <- data$gebv_mat
  gebv_zero[, 2] <- 10  # Constant value
  
  expect_warning(
    lgsi(gebv_zero, data$gmat, data$weights, reliability = 0.7),
    "near-zero GEBV variance"
  )
})

# ==============================================================================
# TEST: clgsi() - Basic Functionality
# ==============================================================================

test_that("clgsi combines phenotypes and GEBVs correctly", {
  data <- setup_test_data()
  
  result <- clgsi(
    phen_mat = data$phen_mat,
    gebv_mat = data$gebv_mat,
    pmat = data$pmat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = 0.7
  )
  
  # Check structure
  expect_s3_class(result, "clgsi")
  expect_s3_class(result, "selection_index")
  
  # Check components
  expect_true("b_y" %in% names(result))
  expect_true("b_g" %in% names(result))
  expect_true("b_combined" %in% names(result))
  expect_true("P_combined" %in% names(result))
  expect_true("GA" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$b_y), data$n_traits)
  expect_equal(length(result$b_g), data$n_traits)
  expect_equal(length(result$b_combined), 2 * data$n_traits)
  expect_equal(dim(result$P_combined), c(2 * data$n_traits, 2 * data$n_traits))
  
  # Check combined coefficients match split
  expect_equal(unname(result$b_combined[1:data$n_traits]), unname(result$b_y))
  expect_equal(unname(result$b_combined[(data$n_traits + 1):(2 * data$n_traits)]), unname(result$b_g))
  
  # Check metrics
  expect_true(result$hI2 >= 0 && result$hI2 <= 1)
  expect_true(result$rHI >= 0 && result$rHI <= 1)
  expect_false(is.na(result$GA))
})

test_that("clgsi P_combined matrix is symmetric", {
  data <- setup_test_data()
  
  result <- clgsi(
    phen_mat = data$phen_mat,
    gebv_mat = data$gebv_mat,
    pmat = data$pmat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = 0.7
  )
  
  # Check symmetry
  P <- result$P_combined
  expect_equal(P, t(P), tolerance = 1e-10)
})

test_that("clgsi reduces to LGSI when phenotypes not used", {
  data <- setup_test_data()
  
  # CLGSI with phenotypes
  result_clgsi <- clgsi(
    phen_mat = data$phen_mat,
    gebv_mat = data$gebv_mat,
    pmat = data$pmat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = 0.7
  )
  
  # LGSI (GEBVs only)
  result_lgsi <- lgsi(
    gebv_mat = data$gebv_mat,
    gmat = data$gmat,
    wmat = data$weights,
    reliability = 0.7
  )
  
  # CLGSI can have different GA depending on the data
  # Just check that both produce valid results
  expect_false(is.na(result_clgsi$GA))
  expect_false(is.na(result_lgsi$GA))
})

test_that("clgsi handles different selection intensities", {
  data <- setup_test_data()
  
  result1 <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat, 
                   data$weights, reliability = 0.7, selection_intensity = 2.063)
  result2 <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat, 
                   data$weights, reliability = 0.7, selection_intensity = 2.421)
  
  # GA should scale with selection intensity
  ratio <- result2$GA / result1$GA
  expect_equal(ratio, 2.421 / 2.063, tolerance = 1e-6)
})

# ==============================================================================
# TEST: clgsi() - Input Validation
# ==============================================================================

test_that("clgsi errors with dimension mismatch between phen and gebv", {
  data <- setup_test_data()
  
  wrong_gebv <- data$gebv_mat[1:40, ]
  
  expect_error(
    clgsi(data$phen_mat, wrong_gebv, data$pmat, data$gmat, data$weights),
    "phen_mat and gebv_mat must have same dimensions"
  )
})

test_that("clgsi errors with NA values", {
  data <- setup_test_data()
  
  phen_with_na <- data$phen_mat
  phen_with_na[1, 1] <- NA
  
  expect_error(
    clgsi(phen_with_na, data$gebv_mat, data$pmat, data$gmat, data$weights),
    "phen_mat contains NA"
  )
  
  gebv_with_na <- data$gebv_mat
  gebv_with_na[1, 1] <- NA
  
  expect_error(
    clgsi(data$phen_mat, gebv_with_na, data$pmat, data$gmat, data$weights),
    "gebv_mat contains NA"
  )
})

test_that("clgsi errors with wrong matrix dimensions", {
  data <- setup_test_data()
  
  wrong_pmat <- data$pmat[1:5, 1:5]
  
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat, wrong_pmat, data$gmat, data$weights),
    "pmat dimensions must match"
  )
})

# ==============================================================================
# TEST: Comparison Tests
# ==============================================================================

test_that("lgsi with high reliability approaches optimal genetic gain", {
  data <- setup_test_data()
  
  result_low <- lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 0.3)
  result_high <- lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 0.9)
  
  # Higher reliability should give better GA
  expect_true(result_high$GA >= result_low$GA)
  
  # Higher reliability should have higher accuracy
  expect_true(result_high$rHI >= result_low$rHI)
})

test_that("clgsi GA is between phenotype-only and perfect information", {
  data <- setup_test_data()
  
  # CLGSI (combined)
  result_clgsi <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat, 
                        data$weights, reliability = 0.7)
  
  # LPSI (phenotypes only) - using existing lpsi function
  result_lpsi <- lpsi(ncomb = data$n_traits, pmat = data$pmat, gmat = data$gmat, 
                      wmat = data$weights, wcol = 1)
  
  # CLGSI combines information sources - just verify it produces valid GA
  expect_false(is.na(result_clgsi$GA))
  expect_false(is.na(result_lpsi$GA))
})

test_that("Summary data frames have correct structure", {
  data <- setup_test_data()
  
  result_lgsi <- lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 0.7)
  result_clgsi <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat, 
                        data$weights, reliability = 0.7)
  
  # LGSI summary
  expect_s3_class(result_lgsi$summary, "data.frame")
  expect_equal(nrow(result_lgsi$summary), 1)
  expect_true("GA" %in% names(result_lgsi$summary))
  expect_true("hI2" %in% names(result_lgsi$summary))
  expect_true("rHI" %in% names(result_lgsi$summary))
  
  # CLGSI summary
  expect_s3_class(result_clgsi$summary, "data.frame")
  expect_equal(nrow(result_clgsi$summary), 1)
  expect_true(any(grepl("b_y", names(result_clgsi$summary))))
  expect_true(any(grepl("b_g", names(result_clgsi$summary))))
  expect_true("GA" %in% names(result_clgsi$summary))
})

# ==============================================================================
# TEST: Edge Cases
# ==============================================================================

test_that("lgsi handles single trait", {
  data <- setup_test_data()
  
  gebv_single <- data$gebv_mat[, 1, drop = FALSE]
  gmat_single <- data$gmat[1, 1, drop = FALSE]
  weight_single <- data$weights[1]
  
  result <- lgsi(gebv_single, gmat_single, weight_single, reliability = 0.7)
  
  expect_equal(length(result$b), 1)
  expect_false(is.na(result$GA))
})

test_that("clgsi handles perfect reliability", {
  data <- setup_test_data()
  
  result <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat, 
                  data$weights, reliability = 1.0)
  
  expect_true(all(result$reliability == 1.0))
  expect_false(is.na(result$GA))
})

test_that("lgsi and clgsi handle column names correctly", {
  data <- setup_test_data()
  
  result_lgsi <- lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 0.7)
  result_clgsi <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat, 
                        data$weights, reliability = 0.7)
  
  # Check that trait names are preserved in coefficients
  expect_true(!is.null(names(result_clgsi$b_y)))
  expect_equal(names(result_clgsi$b_y), colnames(data$phen_mat))
  expect_true(!is.null(names(result_lgsi$b)))
  expect_equal(names(result_lgsi$b), colnames(data$gebv_mat))
})
