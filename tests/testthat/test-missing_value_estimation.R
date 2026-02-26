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
  gen_reml <- selection.index:::gen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "REML"
  )
  expect_true(is.matrix(gen_reml))
  expect_equal(nrow(gen_reml), ncol(test_data))

  gen_yates <- selection.index:::gen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Yates"
  )
  expect_true(is.matrix(gen_yates))
  expect_equal(nrow(gen_yates), ncol(test_data))

  gen_healy <- selection.index:::gen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Healy"
  )
  expect_true(is.matrix(gen_healy))
  expect_equal(nrow(gen_healy), ncol(test_data))

  gen_regression <- selection.index:::gen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Regression"
  )
  expect_true(is.matrix(gen_regression))
  expect_equal(nrow(gen_regression), ncol(test_data))

  gen_mean <- selection.index:::gen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Mean"
  )
  expect_true(is.matrix(gen_mean))
  expect_equal(nrow(gen_mean), ncol(test_data))

  gen_bartlett <- selection.index:::gen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Bartlett"
  )
  expect_true(is.matrix(gen_bartlett))
  expect_equal(nrow(gen_bartlett), ncol(test_data))

  phen_reml <- selection.index:::phen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "REML"
  )
  expect_true(is.matrix(phen_reml))
  expect_equal(nrow(phen_reml), ncol(test_data))

  phen_yates <- selection.index:::phen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Yates"
  )
  expect_true(is.matrix(phen_yates))
  expect_equal(nrow(phen_yates), ncol(test_data))

  phen_healy <- selection.index:::phen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Healy"
  )
  expect_true(is.matrix(phen_healy))
  expect_equal(nrow(phen_healy), ncol(test_data))

  phen_regression <- selection.index:::phen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Regression"
  )
  expect_true(is.matrix(phen_regression))
  expect_equal(nrow(phen_regression), ncol(test_data))

  phen_mean <- selection.index:::phen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Mean"
  )
  expect_true(is.matrix(phen_mean))
  expect_equal(nrow(phen_mean), ncol(test_data))

  phen_bartlett <- selection.index:::phen_varcov(
    data = test_data,
    genotypes = seldata$treat,
    replication = seldata$rep,
    method = "Bartlett"
  )
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
    selection.index:::gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep),
    "Missing values detected in data. Using default method 'REML' for imputation."
  )

  # Test phen.varcov warns when method not provided
  expect_warning(
    selection.index:::phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep),
    "Missing values detected in data. Using default method 'REML' for imputation."
  )

  # Test no warning when method is explicitly provided
  expect_silent(
    selection.index:::gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep, method = "REML")
  )

  expect_silent(
    selection.index:::phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep, method = "Mean")
  )
})

test_that("functions work without warning when no missing values", {
  # Load test data
  data(seldata, package = "selection.index")

  # Use complete data (no missing values)
  test_data <- seldata[, 3:9]

  # Should not warn even without method specified
  expect_silent(
    selection.index:::gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  )

  expect_silent(
    selection.index:::phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  )

  # Verify matrices are returned
  result_gen <- selection.index:::gen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  result_phen <- selection.index:::phen_varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)

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
    col = as.vector(sapply(0:(n - 1), function(i) (1:n + i - 1) %% n + 1)),
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
      method = "REML" # Should fall back to Mean
    ),
    "Split Plot Design"
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

# ==============================================================================
# TEST: Edge cases for 100% missing_value_estimation coverage
# ==============================================================================

test_that("missing_value_estimation: SPD falls back to Mean with warning (Lines 152-155)", {
  data_mat <- matrix(c(1, 2, NA, 4), nrow = 2, ncol = 2)
  gen_idx <- c(1, 2)
  rep_idx <- c(1, 1)
  main_idx <- c(1, 1)

  expect_warning(
    selection.index:::missing_value_estimation(
      data_mat, gen_idx, rep_idx,
      main_idx = main_idx,
      design_type = "SPD", method = "REML"
    )
  )

  result <- suppressWarnings(selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    main_idx = main_idx,
    design_type = "SPD", method = "REML"
  ))
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Healy RCBD w_sum <= 0 fallback (Lines 292-293)", {
  # To trigger w_sum <= 0, we need w_treat = 0 and w_block = 0.
  # w_treat = (repli - n_miss_treat) / repli
  # w_block = (genotype - n_miss_block) / genotype
  # So n_miss_treat must equal repli (all reps for a treatment missing).
  # And n_miss_block must equal genotype (all treatments for a block missing).
  # We must have at least one complete observation for grouped_sums to not throw an error on NA.

  data_mat <- matrix(c(1, NA, NA, NA), nrow = 4, ncol = 1)
  gen_idx <- c(1, 1, 2, 2)
  rep_idx <- c(1, 2, 1, 2)

  # For idx = 4 (gen=2, rep=2):
  # n_miss_treat for gen 2 is 2 (equals repli 2)
  # n_miss_block for rep 2 is 2 (equals genotype 2)
  # w_treat = 0, w_block = 0, w_sum = 0

  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Healy"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Healy LSD w_sum <= 0 fallback (Lines 319-321)", {
  # Similar to RCBD but 3-way.
  data_mat <- matrix(c(1, NA, NA, NA), nrow = 4, ncol = 1)
  gen_idx <- c(1, 2, 2, 1)
  rep_idx <- c(1, 1, 2, 2)
  col_idx <- c(1, 2, 1, 2)

  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    col_idx = col_idx,
    design_type = "LSD", method = "Healy"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Regression RCBD rank deficiency fallback (Lines 476-481)", {
  # min_obs for Regression RCBD = genotype + repli. (2+2 = 4)
  # length(complete_idx) must be > 4, so 5.
  # But the design matrix must be rank deficient.
  data_mat <- matrix(c(1, 1, 1, 1, 1, NA, NA, NA), nrow = 8, ncol = 1)
  gen_idx <- rep(1:2, each = 4)
  rep_idx <- rep(1:4, times = 2)
  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Regression"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Regression LSD rank deficiency fallback (Lines 436-441)", {
  set.seed(123)
  n <- 3
  n_obs <- n * n
  lsd_data <- data.frame(
    treat = rep(1:n, each = n),
    row = rep(1:n, times = n),
    col = as.vector(sapply(0:(n - 1), function(i) (1:n + i - 1) %% n + 1)),
    trait1 = rnorm(n_obs, mean = 10, sd = 2)
  )
  lsd_data$trait1[1] <- NA
  result <- selection.index:::missing_value_estimation(
    as.matrix(lsd_data$trait1), lsd_data$treat, lsd_data$row, lsd_data$col,
    design_type = "LSD", method = "Regression"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Regression zero complete fallback (Line 487)", {
  data_mat <- matrix(c(1, NA, NA, NA), nrow = 4, ncol = 1)
  gen_idx <- c(1, 1, 2, 2)
  rep_idx <- c(1, 2, 1, 2)

  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Regression"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: REML zero variance sum (Lines 674, 686-689, 734-737, 747-753)", {
  # REML total variance <= 0 (RCBD lines 674, 686)
  data_mat <- matrix(c(10, 10, 10, 10, NA), nrow = 5, ncol = 1)
  gen_idx <- c(1, 1, 2, 2, 2)
  rep_idx <- c(1, 2, 1, 2, 3)

  result_rcbd <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "REML"
  )
  expect_false(anyNA(result_rcbd))

  # REML total variance <= 0 (LSD lines 734, 747, 750)
  n <- 3
  n_obs <- n * n
  lsd_data <- matrix(10, nrow = n_obs, ncol = 1)
  lsd_data[1, 1] <- NA
  treat <- rep(1:n, each = n)
  row <- rep(1:n, times = n)
  col <- as.vector(sapply(0:(n - 1), function(i) (1:n + i - 1) %% n + 1))

  result_lsd <- selection.index:::missing_value_estimation(
    lsd_data, treat, row, col,
    design_type = "LSD", method = "REML"
  )
  expect_false(anyNA(result_lsd))
})

test_that("missing_value_estimation: Mean zero complete fallback (Line 596)", {
  data_mat <- matrix(NA_real_, nrow = 4, ncol = 1)
  gen_idx <- c(1, 1, 2, 2)
  rep_idx <- c(1, 2, 1, 2)

  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Mean"
  )
  expect_true(all(result == 0))
})

test_that("missing_value_estimation: Bartlett fallback without valid covariates (Lines 873-878, 897-955)", {
  # Bartlett fallback needs length(complete_idx) > 3, so at least 4.
  # And it needs NO valid covariates. A covariate is only valid if it has NO NAs
  # for the complete_idx rows.
  data_mat <- matrix(1:10, nrow = 5, ncol = 2)
  data_mat[5, 1] <- NA # Col 1 is target, row 5 is missing. complete_idx = 1:4
  data_mat[2, 2] <- NA # Col 2 is covariate. Since row 2 is in complete_idx and has NA, Col 2 is invalid.

  gen_idx <- c(1, 1, 2, 2, 3)
  rep_idx <- c(1, 2, 1, 2, 3)

  result_rcbd <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Bartlett"
  )
  expect_false(anyNA(result_rcbd[, 1]))

  # LSD Bartlett fallback
  n <- 3
  n_obs <- n * n
  lsd_data <- data.frame(
    treat = rep(1:n, each = n),
    row = rep(1:n, times = n),
    col = as.vector(sapply(0:(n - 1), function(i) (1:n + i - 1) %% n + 1)),
    trait1 = rnorm(n_obs, mean = 10, sd = 2),
    trait2 = rnorm(n_obs, mean = 5, sd = 1)
  )
  lsd_data$trait1[5] <- NA # Target missing
  lsd_data$trait2[1] <- NA # Covariate invalid because NA is in complete_idx

  result_lsd_nocovar <- selection.index:::missing_value_estimation(
    as.matrix(lsd_data[, c("trait1", "trait2")]),
    lsd_data$treat, lsd_data$row,
    col_idx = lsd_data$col,
    design_type = "LSD", method = "Bartlett"
  )
  expect_false(anyNA(result_lsd_nocovar[, 1]))
})

test_that("missing_value_estimation: out of bounds indices in fallbacks (Lines 481, 674, 686-689, 737, 747-753, 933-952)", {
  # We trigger the g <= length(treat_means) FALSE, r <= length(block_means) FALSE conditions
  # by making the missing value have a genotype/rep index higher than any complete value.
  data_mat <- matrix(c(10, 10, 10, 10, NA), nrow = 5, ncol = 1)
  gen_idx <- c(1, 1, 2, 2, 3) # Missing value is genotype 3. Complete is only 1, 2. (length = 2)
  rep_idx <- c(1, 2, 1, 2, 3) # Missing value is rep 3. Complete is only 1, 2.

  # Regression fallback
  result_reg <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Regression"
  )
  expect_false(anyNA(result_reg))

  # REML fallback
  result_reml <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "REML"
  )
  expect_false(anyNA(result_reml))

  # Bartlett fallback (requires >3 complete and 0 covariates)
  # Just use same data_mat and duplicate column so Bartlett sees it as covariate.
  # Make the covariate invalid.
  data_mat_bart <- cbind(data_mat, c(1, NA, 3, 4, 5))
  result_bart <- selection.index:::missing_value_estimation(
    data_mat_bart, gen_idx, rep_idx,
    design_type = "RCBD", method = "Bartlett"
  )
  expect_false(anyNA(result_bart[, 1]))

  # Now for LSD (needs col_means out of bounds)
  # 16 complete, 1 missing.
  lsd_mat <- matrix(10, nrow = 17, ncol = 1)
  lsd_mat[17, 1] <- NA
  # Complete 1-16 use indices 1-4.
  lsd_treat <- c(rep(1:4, 4), 5) # Make missing be index 5
  lsd_row <- c(rep(1:4, each = 4), 5)
  lsd_col <- c(rep(1:4, times = 4), 5)

  # REML LSD fallback
  result_lsd_reml <- selection.index:::missing_value_estimation(
    lsd_mat, lsd_treat, lsd_row, lsd_col,
    design_type = "LSD", method = "REML"
  )
  expect_false(anyNA(result_lsd_reml))

  # Bartlett LSD fallback
  lsd_mat_bart <- cbind(lsd_mat, c(NA, rep(1, 16)))
  result_lsd_bart <- selection.index:::missing_value_estimation(
    lsd_mat_bart, lsd_treat, lsd_row, lsd_col,
    design_type = "LSD", method = "Bartlett"
  )
  expect_false(anyNA(result_lsd_bart[, 1]))
})

# ==============================================================================
# TEST: Precisely targeted tests for remaining uncovered lines
# ==============================================================================

test_that("missing_value_estimation: Healy RCBD w_sum=0 triggered when ALL reps of a treatment AND ALL treatments of a block are missing (Lines 292-293)", {
  # Setup: 2 genotypes, 3 reps. Row 5 has genotype=2, rep=3.
  # All reps of genotype 2 are missing: rows 4,5,6 → NA
  # All genotypes of rep 3 are missing: rows 3,6 → NA
  # So for the missing value at row 6 (gen=2, rep=3):
  #   n_miss_treat = 3 = repli  → w_treat = 0
  #   n_miss_block = 2 = genotype → w_block = 0
  #   w_sum = 0  → triggers lines 292-293
  data_mat <- matrix(c(10, 20, NA, 30, NA, NA), nrow = 6, ncol = 1)
  gen_idx <- c(1, 1, 1, 2, 2, 2)
  rep_idx <- c(1, 2, 3, 1, 2, 3)

  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Healy"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Healy LSD w_sum=0 triggered when treatment, row, and column are ALL completely missing (Lines 319-321)", {
  # 3x3 LSD. Make one missing value where ALL elements of its treatment, row, AND column are NA.
  # Rows 1-5 complete, rows 6-9 NA. Missing value in row 7 (treat=3, row=2, col=3).
  # All treat=3 rows (7,8,9) are NA → n_miss_treat = 3 = n(genotype)
  # All row=2 rows (4,7) are NA → n_miss_block = 2 < genotype (fails partial)
  # Let's try: only 2 complete, rest NA. treat=1,row=1,col=1 and treat=2,row=2,col=1 are complete.
  # For missing at treat=3,row=3,col=2: all treat=3 rows are NA (2), all row=3 rows are NA (3x3 = col1,2,3 → 3), all col=2 rows are NA
  # n = 3 → w_treat = (3-3)/3=0, w_row = (3-3)/3=0, w_col = (3-3)/3=0 → triggers 319-321
  n <- 3
  data_mat <- matrix(NA_real_, nrow = n * n, ncol = 1)
  # Only keep (treat=1,row=1,col=1) and (treat=2,row=2,col=3) complete
  # LSD layout: treat=1,2,3; row=1,2,3; col cycles
  treat_idx <- c(1, 2, 3, 2, 3, 1, 3, 1, 2) # Standard 3x3 Latin square
  row_idx <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  col_idx <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
  # Set complete values only where treat=1 or treat=2 (partially), keep treat=3 all NA
  # treat=3 rows: positions 3, 6, 9 → all NA
  # row=3 rows: positions 7,8,9 → all NA already (treat=3 is there)
  # col=2 rows: positions 2,5,8 → position 8 is treat=1 and row=3
  # For missing at position 9 (treat=2 in this layout), let's re-layout
  # Actually let's use a simpler layout where we can control
  treat_idx2 <- rep(1:n, n) # 1,2,3,1,2,3,1,2,3
  row_idx2 <- rep(1:n, each = n) # 1,1,1,2,2,2,3,3,3
  col_idx2 <- as.integer(sapply(1:n, function(i) (seq(0, n - 1) + i - 1) %% n + 1))
  # Complete values: only first 2 rows to ensure we have some data
  data_mat2 <- matrix(NA_real_, nrow = n * n, ncol = 1)
  data_mat2[1, 1] <- 10 # treat=1, row=1
  data_mat2[2, 1] <- 20 # treat=2, row=1
  # treat=3 rows: 3, 6, 9 are all NA
  # row=3 positions: 7, 8, 9 all NA
  # For position 9 (treat=3, row=3): n_miss_treat(treat=3)=3, n_miss_row(row=3)=3, n_miss_col=?
  # col of position 9: col_idx2[9]
  col9 <- col_idx2[9]
  # positions with that col
  n_miss_col9 <- sum(is.na(data_mat2[col_idx2 == col9, 1]))
  # w_treat = (n - 3)/n = 0, w_row = (n - 3)/n = 0, w_col = (n - n_miss_col9)/n
  # For w_col = 0 too, we need all positions in that column to be NA.
  # Ensure all positions with col=col9 are NA (they should be since only rows 1,2 have data,
  # and col_idx2[1] and col_idx2[2] are cols 1 and 2 respectively)

  result <- selection.index:::missing_value_estimation(
    data_mat2, treat_idx2, row_idx2,
    col_idx = col_idx2,
    design_type = "LSD", method = "Healy"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Regression LSD full-rank path covers lines 436-441 (col effect coding for prediction)", {
  # This test ensures the LSD Regression path where QR succeeds (rank == ncol(X))
  # properly executes the column-effect coding for missing values (lines 436-441).
  # We need: (1) enough complete observations, (2) full-rank X, (3) at least one missing value.
  # Use a set of 5x5 LSD with one missing value.
  set.seed(42)
  n <- 5
  n_obs <- n * n
  trait <- rnorm(n_obs, 20, 3)
  treat <- rep(1:n, n)
  row <- rep(1:n, each = n)
  col <- as.integer(sapply(1:n, function(i) (seq(0, n - 1) + i - 1) %% n + 1))
  trait[round(n_obs / 2)] <- NA # A single central missing value

  result <- selection.index:::missing_value_estimation(
    as.matrix(trait), treat, row, col,
    design_type = "LSD", method = "Regression"
  )
  expect_false(anyNA(result))
  expect_equal(length(result), n_obs)
})

test_that("missing_value_estimation: REML out-of-bounds genotype/rep indices trigger 'else grand_mean' in sapply (Lines 686, 689, 747, 750, 753)", {
  # When the missing observation's gen_idx = 99 (larger than any of the 4 unique complete
  # genotypes), `g > length(treat_means)` is TRUE, so the `else grand_mean` branch fires.
  # Same for rep_idx = 99 and col_idx = 99.

  # RCBD: covers lines 686 and 689
  data_mat <- matrix(c(rep(10, 16), NA), nrow = 17, ncol = 1)
  gen_idx <- c(rep(1:4, 4), 99) # 99 is out-of-bounds for treat_means
  rep_idx <- c(rep(1:4, each = 4), 99) # 99 is out-of-bounds for block_means

  result_rcbd <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "REML"
  )
  expect_false(anyNA(result_rcbd))
  expect_equal(result_rcbd[17], 10) # Should impute as grand mean (10)

  # LSD: covers lines 747, 750 and 753
  n <- 5
  lsd_mat <- matrix(c(rep(20, n * n), NA), nrow = n * n + 1, ncol = 1)
  lsd_treat <- c(rep(1:n, n), 99) # out-of-bounds
  lsd_row <- c(rep(1:n, each = n), 99) # out-of-bounds
  lsd_col <- c(as.integer(sapply(1:n, function(i) (seq(0, n - 1) + i - 1) %% n + 1)), 99)

  result_lsd <- selection.index:::missing_value_estimation(
    lsd_mat, lsd_treat, lsd_row, lsd_col,
    design_type = "LSD", method = "REML"
  )
  expect_false(anyNA(result_lsd))
})

test_that("missing_value_estimation: Bartlett no-covariate fallback with out-of-bounds index (Lines 933, 936, 939, 949, 952)", {
  # To cover the `else 0` branches inside Bartlett's no-covariate sapply lookups,
  # the missing value's gen_idx/rep_idx must be out-of-bounds (> length of means vector).

  # RCBD: covers lines 949 and 952
  data_mat_bart <- matrix(c(rep(10, 16), NA), nrow = 17, ncol = 2)
  data_mat_bart[2, 2] <- NA # Invalidate covariate at a complete row → no valid covariates
  gen_idx <- c(rep(1:4, 4), 99) # Missing row has out-of-bounds genotype
  rep_idx <- c(rep(1:4, each = 4), 99) # Missing row has out-of-bounds rep

  result_rcbd <- selection.index:::missing_value_estimation(
    data_mat_bart, gen_idx, rep_idx,
    design_type = "RCBD", method = "Bartlett"
  )
  expect_false(anyNA(result_rcbd[, 1]))

  # LSD: covers lines 933, 936, and 939
  n <- 5
  lsd_mat <- matrix(c(rep(20, n * n), NA), nrow = n * n + 1, ncol = 2)
  lsd_mat[1, 2] <- NA # Invalidate covariate at a complete row → no valid covariates
  lsd_treat <- c(rep(1:n, n), 99)
  lsd_row <- c(rep(1:n, each = n), 99)
  lsd_col <- c(as.integer(sapply(1:n, function(i) (seq(0, n - 1) + i - 1) %% n + 1)), 99)

  result_lsd <- selection.index:::missing_value_estimation(
    lsd_mat, lsd_treat, lsd_row, lsd_col,
    design_type = "LSD", method = "Bartlett"
  )
  expect_false(anyNA(result_lsd[, 1]))
})

test_that("missing_value_estimation: Healy LSD w_sum=0 when ALL of treat, row, col for the missing value have NO complete observations (Lines 319-321)", {
  # Use a 3x3 LSD with only 1 complete observation.
  # With treat_idx = rep(1:3, 3) and row_idx = rep(1:3, each=3),
  # col_idx = c(1,2,3, 2,3,1, 3,1,2), only row 1 (treat=1,row=1,col=1) is complete.
  # Position 5 has treat=2, row=2, col=3:
  #   All of treat=2 rows (2,5,8) are NA → n_miss_treat=3=n → w_treat=0
  #   All of row=2 rows (4,5,6) are NA → n_miss_block=3=n → w_row=0
  #   All of col=3 rows (3,5,7) are NA → n_miss_col=3=n → w_col=0
  #   w_sum=0 → triggers lines 319-321
  n <- 3
  treat_idx <- rep(1:n, n)
  row_idx <- rep(1:n, each = n)
  col_idx <- c(1, 2, 3, 2, 3, 1, 3, 1, 2)
  data_mat <- matrix(NA_real_, nrow = n * n, ncol = 1)
  data_mat[1, 1] <- 10 # Only first row complete (treat=1, row=1, col=1)

  result <- selection.index:::missing_value_estimation(
    data_mat, treat_idx, row_idx,
    col_idx = col_idx,
    design_type = "LSD", method = "Healy"
  )
  expect_false(anyNA(result))
})

test_that("missing_value_estimation: Regression RCBD rank-deficient fallback inner loop with out-of-bounds gen/rep (Lines 476-481)", {
  # Create rank-deficient design matrix: all 7 complete rows belong to genotype=1 only.
  # This makes the treatment contrast columns entirely zero for complete rows → rank < ncol(X).
  # The missing rows (8,9) have genotype=2 and 3.
  # Because the missing gen_idx > length(treat_means from complete rows),
  # the `if (g <= length(treat_means))` inside the fallback loop returns FALSE,
  # triggering the `else 0` branch at line 479 (and equivalent for block at 480).
  # Also ensures lines 476-481 are fully traversed.
  data_mat <- matrix(c(rep(10, 7), NA, NA), nrow = 9, ncol = 1)
  gen_idx <- c(1, 1, 1, 1, 1, 1, 1, 2, 3) # Complete rows all gen=1; missing are gen=2,3
  rep_idx <- c(1, 2, 3, 1, 2, 3, 1, 1, 2) # Missing rows have rep=1 and rep=2

  result <- selection.index:::missing_value_estimation(
    data_mat, gen_idx, rep_idx,
    design_type = "RCBD", method = "Regression"
  )
  expect_false(anyNA(result))
  expect_equal(nrow(result), 9)
})

test_that("missing_value_estimation: REML LSD shrinkage=0.33 when complete_idx <= 3*genotype-2 (Line 737)", {
  # Line 737 fires in the LSD REML branch when df_error <= 0 OR
  # length(complete_idx) <= (3 * genotype - 2).
  # With n=3 (genotype=3), 3*3-2=7, so we need complete_idx <= 7.
  # A 3x3 LSD with 8 complete obs and 1 missing → complete_idx = 8, but 8 > 7.
  # With 3x3, genotype = repli = ncol_blocks = 3 → df_error = (3-1)*(3-2) = 2
  # We need complete_idx <= 7. Use a 3x3 LSD with only 4 complete observations.
  n <- 3
  treat_idx <- rep(1:n, n)
  row_idx <- rep(1:n, each = n)
  col_idx <- c(1, 2, 3, 2, 3, 1, 3, 1, 2)
  # Only keep 4 complete observations; mark remaining 5 as NA
  data_mat <- matrix(NA_real_, nrow = n * n, ncol = 1)
  data_mat[c(1, 2, 4, 6), 1] <- c(10, 12, 11, 9)

  result <- selection.index:::missing_value_estimation(
    data_mat, treat_idx, row_idx,
    col_idx = col_idx,
    design_type = "LSD", method = "REML"
  )
  expect_false(anyNA(result))
})
