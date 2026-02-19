# ==============================================================================
# Tests for Validation Helper Functions (R/helpers-validation.R)
# ==============================================================================

# Note: These are internal validation functions used throughout the package

# ==============================================================================
# TEST: validate_design_args()
# ==============================================================================

test_that("validate_design_args accepts valid RCBD code", {
  result <- selection.index:::validate_design_args(design_type = 1)
  expect_equal(result, 1)
})

test_that("validate_design_args accepts valid LSD code", {
  col_idx <- rep(1:5, times = 5)
  result <- selection.index:::validate_design_args(design_type = 2, col_idx = col_idx)
  expect_equal(result, 2)
})

test_that("validate_design_args accepts valid SPD code", {
  main_idx <- rep(1:4, each = 6)
  result <- selection.index:::validate_design_args(design_type = 3, main_idx = main_idx)
  expect_equal(result, 3)
})

test_that("validate_design_args converts character to code with allow_char", {
  result_rcbd <- selection.index:::validate_design_args("RCBD", allow_char = TRUE)
  expect_equal(result_rcbd, 1)
  
  result_lsd <- selection.index:::validate_design_args("LSD", 
                                                        col_idx = rep(1:5, 5), 
                                                        allow_char = TRUE)
  expect_equal(result_lsd, 2)
  
  result_spd <- selection.index:::validate_design_args("SPD", 
                                                        main_idx = rep(1:3, each = 8),
                                                        allow_char = TRUE)
  expect_equal(result_spd, 3)
})

test_that("validate_design_args errors with invalid design code", {
  expect_error(
    selection.index:::validate_design_args(design_type = 99),
    "must be 1.*2.*3"
  )
  
  expect_error(
    selection.index:::validate_design_args(design_type = 0),
    "must be 1.*2.*3"
  )
})

test_that("validate_design_args errors with invalid character design", {
  expect_error(
    selection.index:::validate_design_args("INVALID", allow_char = TRUE),
    "RCBD.*LSD.*SPD"
  )
})

test_that("validate_design_args errors when LSD missing col_idx", {
  expect_error(
    selection.index:::validate_design_args(design_type = 2),
    "Latin Square.*col_idx"
  )
})

test_that("validate_design_args errors when SPD missing main_idx", {
  expect_error(
    selection.index:::validate_design_args(design_type = 3),
    "Split Plot.*main_idx"
  )
})

test_that("validate_design_args allows RCBD without col_idx or main_idx", {
  # Should not error
  result <- selection.index:::validate_design_args(design_type = 1)
  expect_equal(result, 1)
})

# ==============================================================================
# TEST: validate_indices()
# ==============================================================================

test_that("validate_indices accepts valid indices", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 10)
  
  # Should not error
  expect_silent(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx)
  )
})

test_that("validate_indices errors when genotype length mismatch", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 2)  # Only 20 values
  rep_idx <- rep(1:3, times = 10)
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx),
    "Length of 'genotypes'"
  )
})

test_that("validate_indices errors when replication length mismatch", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 5)  # Only 15 values
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx),
    "Length of 'replications'"
  )
})

test_that("validate_indices errors when columns length mismatch", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 10)
  col_idx <- rep(1:5, times = 5)  # Only 25 values
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx, col_idx = col_idx),
    "Length of 'columns'"
  )
})

test_that("validate_indices errors when main_plots length mismatch", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 10)
  main_idx <- rep(1:4, each = 6)  # Only 24 values
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx, main_idx = main_idx),
    "Length of 'main_plots'"
  )
})

test_that("validate_indices errors with NA in genotypes", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  gen_idx[5] <- NA
  rep_idx <- rep(1:3, times = 10)
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx),
    "genotypes.*NA"
  )
})

test_that("validate_indices errors with NA in replications", {
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 10)
  rep_idx[10] <- NA
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx),
    "replications.*NA"
  )
})

test_that("validate_indices errors with NA in columns", {
  n_obs <- 25
  gen_idx <- rep(1:5, each = 5)
  rep_idx <- rep(1:5, times = 5)
  col_idx <- rep(1:5, times = 5)
  col_idx[12] <- NA
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx, col_idx = col_idx),
    "columns.*NA"
  )
})

test_that("validate_indices errors with too few genotype levels", {
  n_obs <- 10
  gen_idx <- rep(1, times = 10)  # Only 1 level
  rep_idx <- rep(1:2, times = 5)
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx),
    "genotypes.*must have at least 2 unique levels"
  )
})

test_that("validate_indices errors with too few replication levels", {
  n_obs <- 10
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1, times = 10)  # Only 1 level
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx),
    "replications.*must have at least 2 unique levels"
  )
})

test_that("validate_indices errors with too few column levels", {
  n_obs <- 25
  gen_idx <- rep(1:5, each = 5)
  rep_idx <- rep(1:5, times = 5)
  col_idx <- rep(1, times = 25)  # Only 1 level
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx, col_idx = col_idx),
    "columns.*must have at least 2 unique levels"
  )
})

test_that("validate_indices uses custom data_name in errors", {
  n_obs <- 15
  gen_idx <- rep(1:5, each = 2)  # Wrong length
  rep_idx <- rep(1:3, times = 5)
  
  expect_error(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx, data_name = "my_data"),
    "my_data"
  )
})

# ==============================================================================
# TEST: warn_pairwise_psd()
# ==============================================================================

test_that("warn_pairwise_psd returns TRUE for positive definite matrix", {
  # Create a positive definite matrix
  mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  
  result <- selection.index:::warn_pairwise_psd(mat)
  
  expect_true(result)
})

test_that("warn_pairwise_psd returns TRUE for positive semi-definite matrix", {
  # Create a PSD matrix with one zero eigenvalue
  mat <- matrix(c(1, 1, 1, 1), nrow = 2)
  
  # Should return TRUE (within tolerance)
  result <- suppressWarnings(
    selection.index:::warn_pairwise_psd(mat, check_symmetry = FALSE)
  )
  
  expect_true(result)
})

test_that("warn_pairwise_psd warns for non-PSD matrix", {
  # Create a matrix with negative eigenvalue
  mat <- matrix(c(1, 3, 3, 1), nrow = 2)
  
  expect_warning(
    result <- selection.index:::warn_pairwise_psd(mat),
    "not positive semi-definite"
  )
  
  expect_false(result)
})

test_that("warn_pairwise_psd shows eigenvalue information in warning", {
  mat <- matrix(c(1, 3, 3, 1), nrow = 2)
  
  expect_warning(
    selection.index:::warn_pairwise_psd(mat),
    "Minimum eigenvalue"
  )
  
  expect_warning(
    selection.index:::warn_pairwise_psd(mat),
    "Maximum eigenvalue"
  )
})

test_that("warn_pairwise_psd warns for non-square matrix", {
  mat <- matrix(1:6, nrow = 2, ncol = 3)
  
  expect_warning(
    result <- selection.index:::warn_pairwise_psd(mat),
    "not square"
  )
  
  expect_false(result)
})

test_that("warn_pairwise_psd warns for non-symmetric matrix", {
  mat <- matrix(c(4, 1, 2, 3), nrow = 2)  # Not symmetric
  
  expect_warning(
    result <- selection.index:::warn_pairwise_psd(mat, check_symmetry = TRUE),
    "not symmetric"
  )
  
  expect_false(result)
})

test_that("warn_pairwise_psd uses custom matrix name in warnings", {
  mat <- matrix(c(1, 3, 3, 1), nrow = 2)
  
  expect_warning(
    selection.index:::warn_pairwise_psd(mat, mat_name = "TestMatrix"),
    "TestMatrix"
  )
})

test_that("warn_pairwise_psd respects tolerance parameter", {
  # Create a matrix with small negative eigenvalue (~ -0.047)
  mat <- matrix(c(1.0, 0.9, 0.9,
                  0.9, 1.0, 0.5,
                  0.9, 0.5, 1.0), nrow = 3)
  
  # With strict tolerance, should warn
  expect_warning(
    result1 <- selection.index:::warn_pairwise_psd(mat, tolerance = 1e-10),
    "not positive semi-definite"
  )
  
  # With loose tolerance, should pass
  result2 <- selection.index:::warn_pairwise_psd(mat, tolerance = 0.05)
  expect_true(result2)
})

test_that("warn_pairwise_psd handles eigenvalue computation errors", {
  # Create a matrix that might cause issues (e.g., with NaN)
  mat <- matrix(c(1, NaN, NaN, 1), nrow = 2)
  
  expect_warning(
    result <- selection.index:::warn_pairwise_psd(mat),
    "eigenvalue computation failed|not symmetric"
  )
})

# ==============================================================================
# TEST: is_symmetric()
# ==============================================================================

test_that("is_symmetric returns TRUE for symmetric matrix", {
  mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  
  result <- selection.index:::is_symmetric(mat)
  
  expect_true(result)
})

test_that("is_symmetric returns FALSE for non-symmetric matrix", {
  mat <- matrix(c(4, 1, 2, 3), nrow = 2)
  
  result <- selection.index:::is_symmetric(mat)
  
  expect_false(result)
})

test_that("is_symmetric handles tolerance parameter", {
  # Nearly symmetric matrix
  mat <- matrix(c(4, 2.0001, 2, 3), nrow = 2)
  
  # Strict tolerance - not symmetric
  result1 <- selection.index:::is_symmetric(mat, tolerance = 1e-10)
  expect_false(result1)
  
  # Loose tolerance - symmetric
  result2 <- selection.index:::is_symmetric(mat, tolerance = 1e-2)
  expect_true(result2)
})

test_that("is_symmetric works with large matrices", {
  set.seed(123)
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  mat <- (mat + t(mat)) / 2  # Make symmetric
  
  result <- selection.index:::is_symmetric(mat)
  
  expect_true(result)
})

test_that("is_symmetric handles named matrices", {
  mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  dimnames(mat) <- list(c("A", "B"), c("A", "B"))
  
  result <- selection.index:::is_symmetric(mat)
  
  expect_true(result)
})

# ==============================================================================
# TEST: is_zero()
# ==============================================================================

test_that("is_zero returns TRUE for zero", {
  result <- selection.index:::is_zero(0)
  expect_true(result)
})

test_that("is_zero returns TRUE for near-zero within tolerance", {
  result1 <- selection.index:::is_zero(1e-12)
  expect_true(result1)
  
  result2 <- selection.index:::is_zero(-1e-12)
  expect_true(result2)
})

test_that("is_zero returns FALSE for non-zero", {
  result <- selection.index:::is_zero(0.01)
  expect_false(result)
})

test_that("is_zero respects tolerance parameter", {
  value <- 0.001
  
  # With default tolerance (1e-10), should be FALSE
  result1 <- selection.index:::is_zero(value)
  expect_false(result1)
  
  # With loose tolerance (0.01), should be TRUE
  result2 <- selection.index:::is_zero(value, tolerance = 0.01)
  expect_true(result2)
})

test_that("is_zero works with vectors", {
  values <- c(0, 1e-12, 0.01, -1e-12, 1)
  
  results <- sapply(values, selection.index:::is_zero)
  
  expect_equal(results, c(TRUE, TRUE, FALSE, TRUE, FALSE))
})

test_that("is_zero handles negative values correctly", {
  result1 <- selection.index:::is_zero(-0.001, tolerance = 0.01)
  expect_true(result1)
  
  result2 <- selection.index:::is_zero(-0.1, tolerance = 0.01)
  expect_false(result2)
})

test_that("is_zero handles NA and Inf", {
  expect_false(selection.index:::is_zero(NA))
  expect_false(selection.index:::is_zero(Inf))
  expect_false(selection.index:::is_zero(-Inf))
})

# ==============================================================================
# TEST: Integration tests for validation helpers
# ==============================================================================

test_that("validation helpers work together in typical workflow", {
  # Simulate typical validation workflow
  n_obs <- 30
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 10)
  
  # Validate design
  design <- selection.index:::validate_design_args(1)
  expect_equal(design, 1)
  
  # Validate indices
  expect_silent(
    selection.index:::validate_indices(n_obs, gen_idx, rep_idx)
  )
  
  # Check a covariance matrix
  set.seed(999)
  cov_mat <- matrix(rnorm(16), nrow = 4, ncol = 4)
  cov_mat <- (cov_mat + t(cov_mat)) / 2  # Make symmetric
  
  is_sym <- selection.index:::is_symmetric(cov_mat)
  expect_true(is_sym)
  
  # All validations pass - workflow continues
  expect_true(TRUE)
})
