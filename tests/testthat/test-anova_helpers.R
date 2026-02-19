# ==============================================================================
# Tests for ANOVA Helpers (R/anova_helpers.R)
# ==============================================================================

# Note: .calculate_anova() is a deprecated internal function that wraps
# the design_stats_api(). These tests ensure backward compatibility.

test_that(".calculate_anova works for RCBD design", {
  # Create test RCBD data
  set.seed(123)
  n_gen <- 10
  n_rep <- 3
  n_traits <- 4
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  # Suppress deprecation warning for testing
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L  # RCBD
    )
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("GMS" %in% names(result))
  expect_true("EMS" %in% names(result))
  expect_true("MSG" %in% names(result))
  expect_true("MSE" %in% names(result))
  expect_true("DFG" %in% names(result))
  expect_true("DFE" %in% names(result))
  expect_true("n_rep" %in% names(result))
  expect_true("n_gen" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$GMS), n_traits)
  expect_equal(length(result$EMS), n_traits)
  expect_equal(dim(result$MSG), c(n_traits, n_traits))
  expect_equal(dim(result$MSE), c(n_traits, n_traits))
  
  # Check degrees of freedom
  expect_equal(result$DFG, n_gen - 1)
  expect_equal(result$DFE, (n_gen - 1) * (n_rep - 1))
  expect_equal(result$n_gen, n_gen)
  expect_equal(result$n_rep, n_rep)
  
  # Check values are finite
  expect_true(all(is.finite(result$GMS)))
  expect_true(all(is.finite(result$EMS)))
  expect_true(all(is.finite(result$MSG)))
  expect_true(all(is.finite(result$MSE)))
})

test_that(".calculate_anova works for LSD design", {
  # Create test Latin Square data
  set.seed(456)
  n_gen <- 5
  n_rep <- 5
  n_col <- 5
  n_traits <- 3
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  col_idx <- rep(1:n_col, length.out = n_obs)
  
  # Suppress deprecation warning
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      col_idx = col_idx,
      design_type = 2L  # LSD
    )
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("GMS" %in% names(result))
  expect_true("EMS" %in% names(result))
  expect_true("MSG" %in% names(result))
  expect_true("MSE" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$GMS), n_traits)
  expect_equal(length(result$EMS), n_traits)
  expect_equal(dim(result$MSG), c(n_traits, n_traits))
  expect_equal(dim(result$MSE), c(n_traits, n_traits))
  
  # Check degrees of freedom for LSD
  expect_equal(result$DFG, n_gen - 1)
  df_error_lsd <- (n_gen - 1) * (n_rep - 1) - (n_col - 1)
  expect_equal(result$DFE, df_error_lsd)
  
  # Check values are finite
  expect_true(all(is.finite(result$GMS)))
  expect_true(all(is.finite(result$EMS)))
})

test_that(".calculate_anova works for SPD design", {
  # Create test Split Plot data
  set.seed(789)
  n_gen <- 8
  n_rep <- 3
  n_main <- 4
  n_traits <- 3
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  main_idx <- rep(1:n_main, each = (n_obs / n_main))
  
  # Suppress deprecation warning
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      main_idx = main_idx,
      design_type = 3L  # SPD
    )
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true("GMS" %in% names(result))
  expect_true("EMS" %in% names(result))
  expect_true("EMS_MAIN" %in% names(result))
  expect_true("DFE_MAIN" %in% names(result))
  expect_true("n_main" %in% names(result))
  
  # Check dimensions
  expect_equal(length(result$GMS), n_traits)
  expect_equal(length(result$EMS), n_traits)
  expect_equal(length(result$EMS_MAIN), n_traits)
  
  # Check SPD-specific components
  expect_false(is.na(result$DFE_MAIN))
  expect_equal(result$n_main, n_main)
  
  # Check values are finite
  expect_true(all(is.finite(result$GMS)))
  expect_true(all(is.finite(result$EMS)))
  expect_true(all(is.finite(result$EMS_MAIN)))
})

test_that(".calculate_anova issues deprecation warning", {
  set.seed(111)
  n_obs <- 30
  n_traits <- 2
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:10, each = 3)
  rep_idx <- rep(1:3, times = 10)
  
  # Expect deprecation warning
  expect_warning(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    ),
    "deprecated"
  )
  
  expect_warning(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    ),
    "design_stats_api"
  )
})

test_that(".calculate_anova returns matrices for MSG and MSE", {
  set.seed(222)
  n_gen <- 5
  n_rep <- 3
  n_traits <- 3
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    )
  )
  
  # MSG and MSE should be square matrices
  expect_true(is.matrix(result$MSG))
  expect_true(is.matrix(result$MSE))
  expect_equal(nrow(result$MSG), ncol(result$MSG))
  expect_equal(nrow(result$MSE), ncol(result$MSE))
  expect_equal(nrow(result$MSG), n_traits)
})

test_that(".calculate_anova diagonal of MSG equals GMS vector", {
  set.seed(333)
  n_gen <- 6
  n_rep <- 4
  n_traits <- 3
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    )
  )
  
  # Diagonal of MSG should equal GMS vector
  expect_equal(diag(result$MSG), result$GMS)
  expect_equal(diag(result$MSE), result$EMS)
})

test_that(".calculate_anova handles single trait", {
  set.seed(444)
  n_gen <- 8
  n_rep <- 3
  n_traits <- 1
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    )
  )
  
  # Check single trait results
  expect_equal(length(result$GMS), 1)
  expect_equal(length(result$EMS), 1)
  expect_equal(dim(result$MSG), c(1, 1))
  expect_equal(dim(result$MSE), c(1, 1))
  
  # Single value should still be finite
  expect_true(is.finite(result$GMS))
  expect_true(is.finite(result$EMS))
})

test_that(".calculate_anova for SPD has NA for non-SPD fields in RCBD/LSD", {
  set.seed(555)
  n_gen <- 5
  n_rep <- 3
  n_traits <- 2
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  # RCBD should have NA for SPD-specific fields
  result_rcbd <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    )
  )
  
  expect_true(all(is.na(result_rcbd$EMS_MAIN)))
  expect_true(is.na(result_rcbd$DFE_MAIN))
  expect_true(is.na(result_rcbd$n_main))
})

test_that(".calculate_anova returns correct count fields", {
  set.seed(666)
  n_gen <- 7
  n_rep <- 4
  n_traits <- 3
  n_obs <- n_gen * n_rep
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    )
  )
  
  expect_equal(result$n_gen, n_gen)
  expect_equal(result$n_rep, n_rep)
})

test_that(".calculate_anova handles unbalanced data gracefully", {
  set.seed(777)
  # Create slightly unbalanced data (some genotypes have fewer reps)
  gen_idx <- c(rep(1:5, each = 3), rep(6:10, each = 2))
  rep_idx <- c(rep(1:3, times = 5), rep(1:2, times = 5))
  n_obs <- length(gen_idx)
  n_traits <- 2
  
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  
  # Should still compute without error
  result <- suppressWarnings(
    selection.index:::.calculate_anova(
      data_mat = data_mat,
      gen_idx = gen_idx,
      rep_idx = rep_idx,
      design_type = 1L
    )
  )
  
  expect_type(result, "list")
  expect_true(all(is.finite(result$GMS)))
  expect_true(all(is.finite(result$EMS)))
})
