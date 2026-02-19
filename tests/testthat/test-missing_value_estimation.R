# Test the modular missing_value_estimation function

# NOTE: Tests for gen_varcov() and phen_varcov() integration require the full
# package to be loaded (via devtools::test() or R CMD check). They will fail
# when running this file individually with test_file().

test_that("missing_value_estimation works with all three methods", {
  # Create sample RCBD data with missing values
  set.seed(123)
  n_gen <- 5
  n_rep <- 3
  n_traits <- 4
  n_obs <- n_gen * n_rep
  
  # Generate test data
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  colnames(data_mat) <- paste0("trait", 1:n_traits)
  
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  # Introduce missing values
  missing_positions <- c(3, 7, 12, 15)
  data_with_missing <- data_mat
  data_with_missing[missing_positions, 1] <- NA
  data_with_missing[c(5, 10), 2] <- NA
  
  # Test REML method
  result_reml <- selection.index:::missing_value_estimation(data_with_missing, gen_idx, rep_idx, method = "REML")
  expect_true(all(is.finite(result_reml)))
  expect_equal(dim(result_reml), dim(data_with_missing))
  expect_false(any(is.na(result_reml)))
  
  # Test Yates method
  result_yates <- selection.index:::missing_value_estimation(data_with_missing, gen_idx, rep_idx, method = "Yates")
  expect_true(all(is.finite(result_yates)))
  expect_equal(dim(result_yates), dim(data_with_missing))
  expect_false(any(is.na(result_yates)))
  
  # Test Healy & Westmacott method
  result_healy <- selection.index:::missing_value_estimation(data_with_missing, gen_idx, rep_idx, method = "Healy")
  expect_true(all(is.finite(result_healy)))
  expect_equal(dim(result_healy), dim(data_with_missing))
  expect_false(any(is.na(result_healy)))
  
  # Test Regression method
  result_regression <- selection.index:::missing_value_estimation(data_with_missing, gen_idx, rep_idx, method = "Regression")
  expect_true(all(is.finite(result_regression)))
  expect_equal(dim(result_regression), dim(data_with_missing))
  expect_false(any(is.na(result_regression)))
  
  # Test Mean substitution method
  result_mean <- selection.index:::missing_value_estimation(data_with_missing, gen_idx, rep_idx, method = "Mean")
  expect_true(all(is.finite(result_mean)))
  expect_equal(dim(result_mean), dim(data_with_missing))
  expect_false(any(is.na(result_mean)))
  
  # Test Bartlett method
  result_bartlett <- selection.index:::missing_value_estimation(data_with_missing, gen_idx, rep_idx, method = "Bartlett")
  expect_true(all(is.finite(result_bartlett)))
  expect_equal(dim(result_bartlett), dim(data_with_missing))
  expect_false(any(is.na(result_bartlett)))
  
  # Test that non-missing values are preserved
  non_missing_idx <- which(is.finite(data_with_missing[, 3]))
  expect_equal(result_reml[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_yates[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_healy[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_regression[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_mean[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_bartlett[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
})

test_that("missing_value_estimation returns data unchanged when no missing values", {
  set.seed(456)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  result <- selection.index:::missing_value_estimation(data_mat, gen_idx, rep_idx)
  expect_equal(result, data_mat)
})

test_that("gen.varcov and phen.varcov use missing_value_estimation correctly", {
  # Load test data
  data(seldata, package = "selection.index")
  
  # Create data with missing values
  test_data <- seldata[, 3:9]
  test_data[3, 1] <- NA
  test_data[7, 2] <- NA
  test_data[15, 3] <- NA
  
  # Test that functions work with all three methods
  gen_reml <- gen_varcov(data = test_data, 
                         genotypes = seldata$treat, 
                         replication = seldata$rep,
                         method = "REML")
  expect_true(is.matrix(gen_reml))
  expect_equal(nrow(gen_reml), ncol(test_data))
  
  gen_yates <- gen_varcov(data = test_data, 
                          genotypes = seldata$treat, 
                          replication = seldata$rep,
                          method = "Yates")
  expect_true(is.matrix(gen_yates))
  expect_equal(nrow(gen_yates), ncol(test_data))
  
  gen_healy <- gen_varcov(data = test_data, 
                          genotypes = seldata$treat, 
                          replication = seldata$rep,
                          method = "Healy")
  expect_true(is.matrix(gen_healy))
  expect_equal(nrow(gen_healy), ncol(test_data))
  
  gen_regression <- gen_varcov(data = test_data, 
                               genotypes = seldata$treat, 
                               replication = seldata$rep,
                               method = "Regression")
  expect_true(is.matrix(gen_regression))
  expect_equal(nrow(gen_regression), ncol(test_data))
  
  gen_mean <- gen_varcov(data = test_data, 
                         genotypes = seldata$treat, 
                         replication = seldata$rep,
                         method = "Mean")
  expect_true(is.matrix(gen_mean))
  expect_equal(nrow(gen_mean), ncol(test_data))
  
  gen_bartlett <- gen_varcov(data = test_data, 
                             genotypes = seldata$treat, 
                             replication = seldata$rep,
                             method = "Bartlett")
  expect_true(is.matrix(gen_bartlett))
  expect_equal(nrow(gen_bartlett), ncol(test_data))
  
  phen_reml <- phen_varcov(data = test_data, 
                           genotypes = seldata$treat, 
                           replication = seldata$rep,
                           method = "REML")
  expect_true(is.matrix(phen_reml))
  expect_equal(nrow(phen_reml), ncol(test_data))
  
  phen_yates <- phen_varcov(data = test_data, 
                            genotypes = seldata$treat, 
                            replication = seldata$rep,
                            method = "Yates")
  expect_true(is.matrix(phen_yates))
  expect_equal(nrow(phen_yates), ncol(test_data))
  
  phen_healy <- phen_varcov(data = test_data, 
                            genotypes = seldata$treat, 
                            replication = seldata$rep,
                            method = "Healy")
  expect_true(is.matrix(phen_healy))
  expect_equal(nrow(phen_healy), ncol(test_data))
  
  phen_regression <- phen_varcov(data = test_data, 
                                 genotypes = seldata$treat, 
                                 replication = seldata$rep,
                                 method = "Regression")
  expect_true(is.matrix(phen_regression))
  expect_equal(nrow(phen_regression), ncol(test_data))
  
  phen_mean <- phen_varcov(data = test_data, 
                           genotypes = seldata$treat, 
                           replication = seldata$rep,
                           method = "Mean")
  expect_true(is.matrix(phen_mean))
  expect_equal(nrow(phen_mean), ncol(test_data))
  
  phen_bartlett <- phen_varcov(data = test_data, 
                               genotypes = seldata$treat, 
                               replication = seldata$rep,
                               method = "Bartlett")
  expect_true(is.matrix(phen_bartlett))
  expect_equal(nrow(phen_bartlett), ncol(test_data))
})

test_that("functions warn when missing values present but method not specified", {
  # Load test data
  data(seldata, package = "selection.index")
  
  # Create data with missing values
  test_data <- seldata[, 3:9]
  test_data[3, 1] <- NA
  test_data[7, 2] <- NA
  
  # Test gen.varcov warns when method not provided
  expect_warning(
    gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep),
    "Missing values detected in data. Using default method 'REML' for imputation."
  )
  
  # Test phen.varcov warns when method not provided
  expect_warning(
    phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep),
    "Missing values detected in data. Using default method 'REML' for imputation."
  )
  
  # Test no warning when method is explicitly provided
  expect_silent(
    gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep, method = "REML")
  )
  
  expect_silent(
    phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep, method = "Mean")
  )
})

test_that("functions work without warning when no missing values", {
  # Load test data
  data(seldata, package = "selection.index")
  
  # Use complete data (no missing values)
  test_data <- seldata[, 3:9]
  
  # Should not warn even without method specified
  expect_silent(
    gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  )
  
  expect_silent(
    phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  )
  
  # Verify matrices are returned
  result_gen <- gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  result_phen <- phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  
  expect_true(is.matrix(result_gen))
  expect_true(is.matrix(result_phen))
  expect_equal(nrow(result_gen), ncol(test_data))
  expect_equal(nrow(result_phen), ncol(test_data))
})

# ==============================================================================
# TEST: Latin Square Design (LSD) support
# ==============================================================================

test_that("missing_value_estimation works with LSD design", {
  # Create LSD test data
  set.seed(123)
  n <- 5
  n_obs <- n * n
  
  lsd_data <- data.frame(
    treat = rep(1:n, each = n),
    row = rep(1:n, times = n),
    col = as.vector(sapply(0:(n-1), function(i) (1:n + i - 1) %% n + 1)),
    trait1 = rnorm(n_obs, mean = 10, sd = 2),
    trait2 = rnorm(n_obs, mean = 15, sd = 3)
  )
  
  # Introduce missing values
  lsd_data$trait1[c(1, 5, 12)] <- NA
  lsd_data$trait2[c(3, 8)] <- NA
  
  imputed <- selection.index:::missing_value_estimation(
    data_mat = as.matrix(lsd_data[, c("trait1", "trait2")]),
    gen_idx = lsd_data$treat,
    rep_idx = lsd_data$row,
    col_idx = lsd_data$col,
    design_type = "LSD",
    method = "REML"
  )
  
  expect_false(anyNA(imputed))
  expect_equal(dim(imputed), c(n_obs, 2))
})

test_that("missing_value_estimation requires columns for LSD", {
  set.seed(456)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  data_mat[1, 1] <- NA
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  expect_error(
    selection.index:::missing_value_estimation(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = "LSD",
      method = "REML"
    ),
    "col_idx"
  )
})

# ==============================================================================
# TEST: Split Plot Design (SPD) support
# ==============================================================================

test_that("missing_value_estimation works with SPD design", {
  # Create SPD test data
  set.seed(789)
  n_main <- 4
  n_sub <- 6
  n_rep <- 3
  n_obs <- n_main * n_sub * n_rep
  
  spd_data <- data.frame(
    mainplot = rep(1:n_main, each = n_sub * n_rep),
    subplot = rep(rep(1:n_sub, each = n_rep), times = n_main),
    block = rep(1:n_rep, times = n_main * n_sub),
    trait1 = rnorm(n_obs, mean = 12, sd = 2.5)
  )
  
  # Introduce missing values
  spd_data$trait1[c(1, 10, 25, 40)] <- NA
  
  # SPD only supports Mean method
  expect_warning(
    imputed <- selection.index:::missing_value_estimation(
      data_mat = as.matrix(spd_data[, "trait1", drop = FALSE]),
      gen_idx = spd_data$subplot,
      rep_idx = spd_data$block,
      main_idx = spd_data$mainplot,
      design_type = "SPD",
      method = "REML"  # Should fall back to Mean
    ),
    "supported for.*only"
  )
  
  expect_false(anyNA(imputed))
  expect_equal(nrow(imputed), n_obs)
})

test_that("missing_value_estimation requires main_plots for SPD", {
  set.seed(999)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  data_mat[1, 1] <- NA
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  expect_error(
    selection.index:::missing_value_estimation(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = "SPD",
      method = "Mean"
    ),
    "main_idx"
  )
})

# ==============================================================================
# TEST: Input validation
# ==============================================================================

test_that("missing_value_estimation errors with invalid design", {
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  expect_error(
    selection.index:::missing_value_estimation(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = "INVALID"
    ),
    "should be one of"
  )
})

test_that("missing_value_estimation errors with invalid method", {
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  expect_error(
    selection.index:::missing_value_estimation(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      method = "INVALID"
    ),
    "should be one of"
  )
})

# ==============================================================================
# TEST: Edge cases and parameter handling
# ==============================================================================

test_that("missing_value_estimation respects tolerance parameter", {
  set.seed(321)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  data_mat[c(1, 5), 1] <- NA
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  # Strict tolerance
  imputed_strict <- selection.index:::missing_value_estimation(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    method = "Yates",
    tolerance = 1e-8
  )
  
  # Loose tolerance
  imputed_loose <- selection.index:::missing_value_estimation(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    method = "Yates",
    tolerance = 1e-4
  )
  
  # Both should complete without error
  expect_false(anyNA(imputed_strict))
  expect_false(anyNA(imputed_loose))
  expect_equal(dim(imputed_strict), dim(imputed_loose))
})

test_that("missing_value_estimation handles numeric genotypes", {
  set.seed(654)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  data_mat[1, 1] <- NA
  
  # Test with numeric indices
  gen_numeric <- rep(1:5, each = 2)
  rep_numeric <- rep(1:2, times = 5)
  
  imputed_numeric <- selection.index:::missing_value_estimation(
    data_mat = data_mat,
    gen_idx = gen_numeric,
    rep_idx = rep_numeric,
    method = "Mean"
  )
  
  expect_false(anyNA(imputed_numeric))
  expect_equal(dim(imputed_numeric), dim(data_mat))
})

test_that("missing_value_estimation handles high proportion of missing values", {
  set.seed(777)
  data_mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  # Introduce many missing values (40% missing in one column)
  missing_idx <- sample(1:10, size = 4)
  data_mat[missing_idx, 1] <- NA
  
  imputed <- selection.index:::missing_value_estimation(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    method = "REML"
  )
  
  # Should still complete
  expect_false(anyNA(imputed))
})
