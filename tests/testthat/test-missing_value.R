# Tests for estimate_missing_values() exported function
# This tests the user-facing wrapper, not the internal implementation

# ==============================================================================
# TEST: Basic functionality - RCBD with all methods
# ==============================================================================

test_that("estimate_missing_values works with RCBD and all methods", {
  data(seldata)

  # Create test data with missing values
  test_data <- seldata[, 3:6]
  test_data[c(1, 10, 25), 1] <- NA
  test_data[c(5, 15), 2] <- NA

  # Test all methods
  methods <- c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")

  for (method in methods) {
    result <- estimate_missing_values(
      data = test_data,
      genotypes = seldata$treat,
      replications = seldata$rep,
      design = "RCBD",
      method = method
    )

    expect_false(anyNA(result), info = paste("Method:", method))
    expect_equal(dim(result), dim(test_data), info = paste("Method:", method))
    expect_true(all(is.finite(result)), info = paste("Method:", method))

    # Check that non-missing values are preserved
    non_missing_idx <- which(is.finite(test_data[, 3]))
    expect_equal(result[non_missing_idx, 3], test_data[non_missing_idx, 3],
      info = paste("Method:", method)
    )
  }
})

test_that("estimate_missing_values works with LSD and all methods", {
  # Create Latin Square Design test data
  set.seed(789)
  n <- 5 # 5x5 Latin square
  lsd_data <- matrix(rnorm(n * n * 3, mean = 20, sd = 5), nrow = n * n, ncol = 3)
  colnames(lsd_data) <- c("trait1", "trait2", "trait3")

  gen_vec <- rep(1:n, times = n)
  row_vec <- rep(1:n, each = n)
  col_vec <- rep(1:n, length.out = n * n)

  # Introduce missing values
  lsd_data[c(2, 7, 15), 1] <- NA
  lsd_data[c(10, 20), 2] <- NA

  methods <- c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")

  for (method in methods) {
    result <- estimate_missing_values(
      data = lsd_data,
      genotypes = gen_vec,
      replications = row_vec,
      columns = col_vec,
      design = "LSD",
      method = method
    )

    expect_false(anyNA(result), info = paste("LSD Method:", method))
    expect_equal(dim(result), dim(lsd_data), info = paste("LSD Method:", method))
    expect_true(all(is.finite(result)), info = paste("LSD Method:", method))
  }
})

test_that("estimate_missing_values works with SPD design", {
  # Create Split Plot Design test data
  set.seed(456)
  n_blocks <- 3
  n_main <- 4
  n_sub <- 5
  n_obs <- n_blocks * n_main * n_sub

  spd_data <- matrix(rnorm(n_obs * 3, mean = 15, sd = 4), nrow = n_obs, ncol = 3)
  colnames(spd_data) <- c("yield", "height", "protein")

  block_vec <- rep(1:n_blocks, each = n_main * n_sub)
  main_vec <- rep(rep(1:n_main, each = n_sub), times = n_blocks)
  sub_vec <- rep(1:n_sub, times = n_blocks * n_main)

  # Introduce missing values
  spd_data[c(5, 15, 25, 35), 1] <- NA
  spd_data[c(10, 30), 2] <- NA

  # SPD only supports Mean method
  result <- estimate_missing_values(
    data = spd_data,
    genotypes = sub_vec,
    replications = block_vec,
    main_plots = main_vec,
    design = "SPD",
    method = "Mean"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(spd_data))
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# TEST: Input type handling
# ==============================================================================

test_that("estimate_missing_values handles data.frame input", {
  data(seldata)

  # Use data frame instead of matrix
  test_df <- seldata[, 3:5]
  test_df[c(1, 10), 1] <- NA

  result <- estimate_missing_values(
    data = test_df,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Yates"
  )

  expect_false(anyNA(result))
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(test_df))
})

test_that("estimate_missing_values handles character genotype labels", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[c(1, 10), 1] <- NA

  # Convert numeric genotypes to character
  gen_char <- paste0("G", seldata$treat)
  rep_char <- paste0("R", seldata$rep)

  result <- estimate_missing_values(
    data = test_data,
    genotypes = gen_char,
    replications = rep_char,
    design = "RCBD",
    method = "Mean"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
})

test_that("estimate_missing_values handles factor genotype labels", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[c(1, 10), 1] <- NA

  # Convert to factors
  gen_factor <- factor(seldata$treat)
  rep_factor <- factor(seldata$rep)

  result <- estimate_missing_values(
    data = test_data,
    genotypes = gen_factor,
    replications = rep_factor,
    design = "RCBD",
    method = "Mean"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
})

# ==============================================================================
# TEST: Column name preservation
# ==============================================================================

test_that("estimate_missing_values preserves column names", {
  data(seldata)

  test_data <- seldata[, 3:5]
  original_names <- colnames(test_data)
  test_data[c(1, 10), 1] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Mean"
  )

  expect_equal(colnames(result), original_names)
})

test_that("estimate_missing_values preserves custom column names", {
  data(seldata)

  test_data <- seldata[, 3:5]
  colnames(test_data) <- c("Yield", "Height", "Protein")
  test_data[c(1, 10), 1] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Mean"
  )

  expect_equal(colnames(result), c("Yield", "Height", "Protein"))
})

test_that("estimate_missing_values handles data without column names", {
  data(seldata)

  test_data <- as.matrix(seldata[, 3:5])
  colnames(test_data) <- NULL # Remove column names
  test_data[c(1, 10), 1] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Mean"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
  expect_null(colnames(result))
})

# ==============================================================================
# TEST: Design-specific validation
# ==============================================================================

test_that("estimate_missing_values requires columns for LSD", {
  set.seed(123)
  lsd_data <- matrix(rnorm(25 * 3), nrow = 25, ncol = 3)
  lsd_data[1, 1] <- NA

  gen_vec <- rep(1:5, times = 5)
  row_vec <- rep(1:5, each = 5)

  expect_error(
    estimate_missing_values(
      data = lsd_data,
      genotypes = gen_vec,
      replications = row_vec,
      design = "LSD",
      method = "Mean"
    ),
    "col_idx"
  )
})

test_that("estimate_missing_values requires main_plots for SPD", {
  set.seed(123)
  spd_data <- matrix(rnorm(60 * 3), nrow = 60, ncol = 3)
  spd_data[1, 1] <- NA

  block_vec <- rep(1:3, each = 20)
  sub_vec <- rep(1:5, times = 12)

  expect_error(
    estimate_missing_values(
      data = spd_data,
      genotypes = sub_vec,
      replications = block_vec,
      design = "SPD",
      method = "Mean"
    ),
    "main_idx"
  )
})

test_that("estimate_missing_values validates design argument", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[1, 1] <- NA

  expect_error(
    estimate_missing_values(
      data = test_data,
      genotypes = seldata$treat,
      replications = seldata$rep,
      design = "InvalidDesign",
      method = "Mean"
    ),
    "should be one of"
  )
})

test_that("estimate_missing_values validates method argument", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[1, 1] <- NA

  expect_error(
    estimate_missing_values(
      data = test_data,
      genotypes = seldata$treat,
      replications = seldata$rep,
      design = "RCBD",
      method = "InvalidMethod"
    ),
    "should be one of"
  )
})

# ==============================================================================
# TEST: Method fallback for SPD
# ==============================================================================

test_that("estimate_missing_values falls back to Mean for SPD with unsupported methods", {
  set.seed(456)
  n_blocks <- 3
  n_main <- 4
  n_sub <- 5
  n_obs <- n_blocks * n_main * n_sub

  spd_data <- matrix(rnorm(n_obs * 3), nrow = n_obs, ncol = 3)
  spd_data[c(5, 15), 1] <- NA

  block_vec <- rep(1:n_blocks, each = n_main * n_sub)
  main_vec <- rep(rep(1:n_main, each = n_sub), times = n_blocks)
  sub_vec <- rep(1:n_sub, times = n_blocks * n_main)

  # Test that unsupported methods issue warning and fall back to Mean
  unsupported_methods <- c("REML", "Yates", "Healy", "Regression", "Bartlett")

  for (method in unsupported_methods) {
    expect_warning(
      result <- estimate_missing_values(
        data = spd_data,
        genotypes = sub_vec,
        replications = block_vec,
        main_plots = main_vec,
        design = "SPD",
        method = method
      ),
      "Mean.*SPD"
    )

    expect_false(anyNA(result), info = paste("SPD fallback for:", method))
  }
})

# ==============================================================================
# TEST: Tolerance parameter
# ==============================================================================

test_that("estimate_missing_values respects tolerance parameter", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[c(1, 10, 25), 1] <- NA
  test_data[c(5, 15), 2] <- NA

  # Test with different tolerance values
  result_strict <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Yates",
    tolerance = 1e-8
  )

  result_loose <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Yates",
    tolerance = 1e-4
  )

  expect_false(anyNA(result_strict))
  expect_false(anyNA(result_loose))

  # Both should produce valid results
  expect_true(all(is.finite(result_strict)))
  expect_true(all(is.finite(result_loose)))
})

# ==============================================================================
# TEST: Edge cases
# ==============================================================================

test_that("estimate_missing_values handles data with no missing values", {
  data(seldata)

  # No missing values
  test_data <- seldata[, 3:5]

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Mean"
  )

  expect_equal(result, as.matrix(test_data))
})

test_that("estimate_missing_values handles data with many missing values", {
  data(seldata)

  test_data <- seldata[, 3:5]
  # Introduce many missing values (about 30%)
  set.seed(999)
  missing_idx <- sample(seq_len(nrow(test_data)), size = floor(0.3 * nrow(test_data)))
  test_data[missing_idx, 1] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "REML"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
})

test_that("estimate_missing_values handles single trait", {
  data(seldata)

  test_data <- seldata[, 3, drop = FALSE]
  test_data[c(1, 10, 25), 1] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Mean"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
})

test_that("estimate_missing_values handles missing values across all traits", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[1, ] <- NA # All traits missing for one observation
  test_data[10, 1] <- NA
  test_data[15, 2] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Mean"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
})

# ==============================================================================
# TEST: Consistency with internal function
# ==============================================================================

test_that("estimate_missing_values produces same results as internal function", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[c(1, 10, 25), 1] <- NA
  test_data[c(5, 15), 2] <- NA

  # Call exported wrapper
  result_wrapper <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "Yates"
  )

  # Call internal function directly
  gen_idx <- as.integer(factor(seldata$treat))
  rep_idx <- as.integer(factor(seldata$rep))

  result_internal <- selection.index:::missing_value_estimation(
    data_mat = as.matrix(test_data),
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    design_type = "RCBD",
    method = "Yates"
  )

  expect_equal(result_wrapper, result_internal)
})

# ==============================================================================
# TEST: Default arguments
# ==============================================================================

test_that("estimate_missing_values uses default design and method", {
  data(seldata)

  test_data <- seldata[, 3:5]
  test_data[c(1, 10), 1] <- NA

  # Test with defaults (RCBD, REML)
  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
})

# ==============================================================================
# TEST: Integration with real data patterns
# ==============================================================================

test_that("estimate_missing_values works with realistic missing patterns", {
  data(seldata)

  test_data <- seldata[, 3:7]

  # Simulate realistic missing patterns:
  # - Random missing values
  # - Entire observation missing
  # - Systematic missing (e.g., measurement failure in one block)

  set.seed(2020)
  # Random missing
  test_data[sample(seq_len(nrow(test_data)), 5), sample(seq_len(ncol(test_data)), 1)] <- NA

  # Entire observation
  test_data[10, ] <- NA

  # Systematic (all observations from block 1, trait 2)
  block1_idx <- which(seldata$rep == 1)
  test_data[block1_idx[1:3], 2] <- NA

  result <- estimate_missing_values(
    data = test_data,
    genotypes = seldata$treat,
    replications = seldata$rep,
    design = "RCBD",
    method = "REML"
  )

  expect_false(anyNA(result))
  expect_equal(dim(result), dim(test_data))
  expect_true(all(is.finite(result)))
})

# ==============================================================================
# TEST: All method combinations for RCBD
# ==============================================================================

test_that("estimate_missing_values works with all RCBD method combinations", {
  data(seldata)

  test_data <- as.matrix(seldata[, 3:5])
  test_data[c(1, 5, 10, 15, 20), 1] <- NA
  test_data[c(2, 8, 14), 2] <- NA

  methods <- c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")

  results <- list()

  for (method in methods) {
    result <- estimate_missing_values(
      data = test_data,
      genotypes = seldata$treat,
      replications = seldata$rep,
      design = "RCBD",
      method = method
    )

    results[[method]] <- result

    # All methods should produce complete results
    expect_false(anyNA(result), info = paste("Method:", method))
    expect_equal(dim(result), dim(test_data), info = paste("Method:", method))

    # Check that non-missing values are identical to original
    non_missing_mask <- is.finite(test_data)
    expect_equal(result[non_missing_mask], test_data[non_missing_mask],
      info = paste("Method:", method)
    )
  }

  # All methods should produce finite results
  for (method in methods) {
    expect_true(all(is.finite(results[[method]])),
      info = paste("Method:", method, "- all finite")
    )
  }
})

# ==============================================================================
# TEST: All method combinations for LSD
# ==============================================================================

test_that("estimate_missing_values works with all LSD method combinations", {
  # Create 6x6 Latin Square Design
  set.seed(555)
  n <- 6
  lsd_data <- matrix(rnorm(n * n * 4, mean = 25, sd = 6), nrow = n * n, ncol = 4)
  colnames(lsd_data) <- paste0("T", 1:4)

  gen_vec <- rep(1:n, times = n)
  row_vec <- rep(1:n, each = n)
  col_vec <- rep(1:n, length.out = n * n)

  # Introduce missing values
  lsd_data[c(2, 10, 20, 30), 1] <- NA
  lsd_data[c(5, 15, 25), 2] <- NA
  lsd_data[c(8, 18), 3] <- NA

  methods <- c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")

  for (method in methods) {
    result <- estimate_missing_values(
      data = lsd_data,
      genotypes = gen_vec,
      replications = row_vec,
      columns = col_vec,
      design = "LSD",
      method = method
    )

    expect_false(anyNA(result), info = paste("LSD Method:", method))
    expect_equal(dim(result), dim(lsd_data), info = paste("LSD Method:", method))

    # Check preservation of non-missing values
    non_missing_mask <- is.finite(lsd_data)
    expect_equal(result[non_missing_mask], lsd_data[non_missing_mask],
      info = paste("LSD Method:", method)
    )
  }
})
