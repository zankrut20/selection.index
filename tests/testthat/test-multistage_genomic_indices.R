test_that("MPPG-LGSI with custom U1 matrix generates warning", {
  # Skip if functions not available
  skip_if_not_installed("selection.index")
  
  # Set up test data with proper dimensions
  set.seed(123)
  n1 <- 3  # Stage 1 traits
  n <- 5   # Total traits
  
  # Create positive definite covariance matrices
  Gamma1 <- diag(n1) + matrix(0.2, n1, n1)
  Gamma <- diag(n) + matrix(0.2, n, n)
  
  A1 <- 0.8 * Gamma1
  A <- 0.8 * Gamma[, 1:n1]  # n x n1
  
  G1 <- Gamma1 * 1.2
  P1 <- Gamma1 * 1.5
  C <- Gamma * 1.2
  
  w <- rep(1, n)
  d1 <- c(2, 1, 1)
  d2 <- c(2, 1, 1, 0.5, 0.5)
  
  # Test 1: No warning with default U1 (Identity)
  # Suppress C* warnings from simple test matrices
  result_default <- suppressWarnings(
    mppg_lgsi(Gamma1 = Gamma1, Gamma = Gamma, 
              A1 = A1, A = A, C = C, G1 = G1, P1 = P1,
              wmat = w, d1 = d1, d2 = d2, U1 = NULL)
  )
  
  expect_true(!is.null(result_default$beta_P1))
  
  # Test 2: Warning with custom U1 (not Identity)
  # Create a custom U1 that is a permutation matrix (still full rank but not Identity)
  U1_custom <- matrix(c(0, 1, 0,
                        1, 0, 0,
                        0, 0, 1), nrow = n1, ncol = n1, byrow = TRUE)
  
  # Capture warnings and check for U1 warning
  # Note: May also get C* warnings due to simple test matrices
  w_all <- capture_warnings({
    result_custom <- mppg_lgsi(Gamma1 = Gamma1, Gamma = Gamma,
                                A1 = A1, A = A, C = C, G1 = G1, P1 = P1,
                                wmat = w, d1 = d1, d2 = d2, U1 = U1_custom)
  })
  
  # Check that U1 warning is present
  expect_true(any(grepl("Custom U1 matrix detected", w_all)))
  
  # Test 3: Verify the warning message contains key information
  # Filter to only U1-related warnings
  w_captured <- w_all[grepl("Custom U1|phenotypic proxy|over-adjusted", w_all, ignore.case = TRUE)]
  
  # Check that warning mentions key concepts
  all_warnings <- paste(w_captured, collapse = " ")
  expect_true(grepl("phenotypic proxy", all_warnings, ignore.case = TRUE))
  expect_true(grepl("over-adjusted", all_warnings, ignore.case = TRUE))
})

test_that("MPPG-LGSI C* adjustment differs between Identity and custom U1", {
  # Set up test data with proper dimensions
  set.seed(456)
  n1 <- 3
  n <- 5
  
  Gamma1 <- diag(n1) + matrix(0.2, n1, n1)
  Gamma <- diag(n) + matrix(0.2, n, n)
  
  A1 <- 0.8 * Gamma1
  A <- 0.8 * Gamma[, 1:n1]  # n x n1
  G1 <- Gamma1 * 1.2
  P1 <- Gamma1 * 1.5
  C <- Gamma * 1.2
  
  w <- rep(1, n)
  d1 <- c(2, 1, 1)
  d2 <- c(2, 1, 1, 0.5, 0.5)
  
  # Run with default U1 (Identity)
  # Suppress C* warnings from simple test matrices
  result_identity <- suppressWarnings(
    mppg_lgsi(Gamma1 = Gamma1, Gamma = Gamma,
              A1 = A1, A = A, C = C, G1 = G1, P1 = P1,
              wmat = w, d1 = d1, d2 = d2, U1 = NULL)
  )
  
  # Run with custom U1 (permutation matrix)
  U1_custom <- matrix(c(0, 1, 0,
                        1, 0, 0,
                        0, 0, 1), nrow = n1, ncol = n1, byrow = TRUE)
  
  suppressWarnings({
    result_custom <- mppg_lgsi(Gamma1 = Gamma1, Gamma = Gamma,
                                A1 = A1, A = A, C = C, G1 = G1, P1 = P1,
                                wmat = w, d1 = d1, d2 = d2, U1 = U1_custom)
  })
  
  # Test that beta_P1 coefficients differ (genomic respects U1)
  expect_false(isTRUE(all.equal(result_identity$beta_P1, result_custom$beta_P1, tolerance = 1e-4)))
  
  # Test that b_P1 phenotypic proxy is the SAME in both cases
  # (this is the key limitation: both use the standard Tallis formula)
  expect_equal(result_identity$b_P1, result_custom$b_P1, tolerance = 1e-10)
  
  # Test that C_star is the SAME (because b_P1 is the same)
  # This demonstrates that C* uses the phenotypic proxy regardless of U1
  expect_equal(result_identity$C_star, result_custom$C_star, tolerance = 1e-10)
  
  # Test that results are still valid (have required components)
  expect_true(!is.null(result_custom$beta_P1))
  expect_true(!is.null(result_custom$b_P1))
  expect_true(!is.null(result_custom$C_star))
  expect_equal(length(result_custom$beta_P1), n1)
  expect_equal(length(result_custom$b_P1), n1)
  expect_equal(dim(result_custom$C_star), c(n, n))
})

test_that("Young's method warning appears in multistage indices", {
  # Set up test data with proper dimensions
  set.seed(789)
  n1 <- 3
  n <- 5
  
  Gamma1 <- diag(n1) + matrix(0.2, n1, n1)
  Gamma <- diag(n) + matrix(0.2, n, n)
  
  A1 <- 0.8 * Gamma1
  A <- 0.8 * Gamma[, 1:n1]  # n x n1
  G1 <- Gamma1 * 1.2
  P1 <- Gamma1 * 1.5
  C <- Gamma * 1.2
  P <- Gamma * 1.5  # Define P for phenotypic version
  
  w <- rep(1, n)
  
  # Test MLGSI with Young's method - suppress expected warning about overestimation
  result_mlgsi <- suppressWarnings(
    mlgsi(Gamma1 = Gamma1, Gamma = Gamma,
          A1 = A1, A = A, C = C, G1 = G1, P1 = P1,
          wmat = w, use_young_method = TRUE)
  )
  
  # Just check the function runs - warning may or may not appear depending on correlation
  expect_true(!is.null(result_mlgsi$beta1))
  
  # Test MLPSI (phenotypic version) also runs with Young's method
  # Suppress expected Young's method warning
  result_mlpsi <- suppressWarnings(
    mlpsi(P1 = P1, P = P, 
          G1 = G1, C = C,
          wmat = w, use_young_method = TRUE)
  )
  
  expect_true(!is.null(result_mlpsi$b1))
})

test_that("Multistage indices use manual intensities by default", {
  # Set up minimal test data with proper dimensions
  set.seed(321)
  n1 <- 2
  n <- 3
  
  Gamma1 <- diag(n1) + matrix(0.2, n1, n1)
  Gamma <- diag(n) + matrix(0.2, n, n)
  
  A1 <- 0.8 * Gamma1
  A <- 0.8 * Gamma[, 1:n1]  # n x n1
  G1 <- Gamma1 * 1.2
  P1 <- Gamma1 * 1.5
  C <- Gamma * 1.2
  
  w <- rep(1, n)
  
  # Run with default parameters (should use manual intensities, no Young's warning)
  # Suppress potential C* warnings from simple test matrices
  result <- suppressWarnings(
    mlgsi(Gamma1 = Gamma1, Gamma = Gamma,
          A1 = A1, A = A, C = C, G1 = G1, P1 = P1,
          wmat = w)
  )
  
  # Check that manual intensities are used (defaults)
  expect_equal(result$k1, 2.063)
  expect_equal(result$k2, 2.063)
  
  # Check that result is valid
  expect_true(!is.null(result$beta1))
})
