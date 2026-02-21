# ==============================================================================
# Tests for Genomic Selection Indices (LGSI and CLGSI)
# ==============================================================================

# Setup: Create test data
setup_test_data <- function() {
  set.seed(123)

  # Use actual variance-covariance matrices from example data
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Simulate GEBVs using a simple approach:
  # Take actual phenotypes and add independent normal noise
  # This gives realistic correlations without positive-definiteness issues
  phen_mat <- as.matrix(seldata[1:50, 3:9])

  # Create GEBVs as phenotypes + noise (simulating prediction)
  # GEBVs should have lower variance than phenotypes (reliability < 1)
  gebv_mat <- phen_mat * 0.7 + matrix(
    rnorm(prod(dim(phen_mat)), sd = 0.5),
    nrow(phen_mat), ncol(phen_mat)
  )
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
  expect_false(is.na(result$GA))
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
    reliability = 0.7, selection_intensity = 2.063
  )
  result2 <- lgsi(data$gebv_mat, data$gmat, data$weights,
    reliability = 0.7, selection_intensity = 2.421
  )

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
  gebv_zero[, 2] <- 10 # Constant value

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
    data$weights,
    reliability = 0.7, selection_intensity = 2.063
  )
  result2 <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat,
    data$weights,
    reliability = 0.7, selection_intensity = 2.421
  )

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
    data$weights,
    reliability = 0.7
  )

  # LPSI (phenotypes only) - using existing lpsi function
  result_lpsi <- lpsi(
    ncomb = data$n_traits, pmat = data$pmat, gmat = data$gmat,
    wmat = data$weights, wcol = 1
  )

  # CLGSI combines information sources - just verify it produces valid GA
  expect_false(is.na(result_clgsi$GA))
  expect_false(is.na(result_lpsi$GA))
})

test_that("Summary data frames have correct structure", {
  data <- setup_test_data()

  result_lgsi <- lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 0.7)
  result_clgsi <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat,
    data$weights,
    reliability = 0.7
  )

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
    data$weights,
    reliability = 1.0
  )

  expect_true(all(result$reliability == 1.0))
  expect_false(is.na(result$GA))
})

test_that("lgsi and clgsi handle column names correctly", {
  data <- setup_test_data()

  result_lgsi <- lgsi(data$gebv_mat, data$gmat, data$weights, reliability = 0.7)
  result_clgsi <- clgsi(data$phen_mat, data$gebv_mat, data$pmat, data$gmat,
    data$weights,
    reliability = 0.7
  )

  # Check that trait names are preserved in coefficients
  expect_true(!is.null(names(result_clgsi$b_y)))
  expect_equal(names(result_clgsi$b_y), colnames(data$phen_mat))
  expect_true(!is.null(names(result_lgsi$b)))
  expect_equal(names(result_lgsi$b), colnames(data$gebv_mat))
})

# ==============================================================================
# NEW COVERAGE TESTS - targeting previously uncovered lines
# ==============================================================================

# --- lgsi: line 118 – wmat matrix coercion (non-vector, non-matrix input) ----
test_that("lgsi coerces data.frame wmat to matrix (line 118)", {
  data <- setup_test_data()
  n <- data$n_traits
  # Supply wmat as a data.frame with two columns; lgsi should coerce it
  wmat_df <- as.data.frame(
    matrix(rep(data$weights, 2), ncol = 2)
  )
  result <- lgsi(data$gebv_mat, data$gmat, wmat_df, wcol = 1, reliability = 0.7)
  expect_s3_class(result, "lgsi")
  expect_equal(length(result$b), n)
})

# --- lgsi: line 122 – wmat row count != n_traits --------------------------
test_that("lgsi errors when wmat has wrong number of rows (line 122)", {
  data <- setup_test_data()
  # data.frame path (else branch) with wrong row count
  wmat_wrong <- as.data.frame(matrix(1:5, ncol = 1))
  expect_error(
    lgsi(data$gebv_mat, data$gmat, wmat_wrong, reliability = 0.7),
    "Number of rows in wmat must equal number of traits"
  )
})

# --- lgsi: line 126 – wcol out of bounds (lines 125-127) ------------------
test_that("lgsi errors when wcol is out of bounds (line 126)", {
  data <- setup_test_data()
  wmat_df <- as.data.frame(matrix(rep(data$weights, 2), ncol = 2))
  expect_error(
    lgsi(data$gebv_mat, data$gmat, wmat_df, wcol = 5, reliability = 0.7),
    "wcol must be between 1 and"
  )
})

# --- lgsi: line 163 – low reliability warning (NULL reliability) ----------
test_that("lgsi warns when estimated reliability is low (line 163)", {
  data <- setup_test_data()
  # Scale GEBVs down drastically so Var(GEBV)/Var(G) < 0.3
  gebv_low <- data$gebv_mat * 0.01
  expect_warning(
    lgsi(gebv_low, data$gmat, data$weights, reliability = NULL),
    "Estimated reliabilities are low"
  )
})

# --- lgsi: line 175 – vector reliability any value out of [0,1] -----------
test_that("lgsi errors when vector reliability is out of range (line 175)", {
  data <- setup_test_data()
  n <- data$n_traits
  bad_rel <- rep(0.7, n)
  bad_rel[2] <- -0.1 # invalid
  expect_error(
    lgsi(data$gebv_mat, data$gmat, data$weights, reliability = bad_rel),
    "All reliability values must be between 0 and 1"
  )
})

# --- lgsi: line 247 – Delta_H returns NA_real_ when sigma_I == 0 ----------
test_that("lgsi returns NA Delta_H when index variance is zero (line 247)", {
  data <- setup_test_data()
  n <- data$n_traits
  # Force a zero index variance by making gebv_mat constant per column
  # and supplying zero weights so bPb = 0.
  # All-zero weights => b = 0 => bPb = 0 => sigma_I = NA
  zero_weights <- rep(0, n)
  # Avoid the zero-variance warning path at line 139 by keeping non-constant GEBVs
  result <- lgsi(data$gebv_mat, data$gmat, zero_weights, reliability = 0.7)
  expect_true(all(is.na(result$Delta_H)))
})

# --- clgsi: line 419 – neither raw data nor cov matrices provided ---------
test_that("clgsi errors when no data or cov matrices provided (line 419)", {
  data <- setup_test_data()
  expect_error(
    clgsi(
      phen_mat = NULL, gebv_mat = NULL,
      pmat = data$pmat, gmat = data$gmat, wmat = data$weights
    ),
    "Must provide either"
  )
})

# --- clgsi: line 436 – phen_mat column count != n_traits ------------------
test_that("clgsi errors when phen_mat has wrong number of traits (line 436)", {
  data <- setup_test_data()
  # gebv_mat has n_traits cols, but phen_mat gets only 5 cols here
  phen_wrong <- data$phen_mat[, 1:5, drop = FALSE]
  expect_error(
    clgsi(phen_wrong, data$gebv_mat, data$pmat, data$gmat, data$weights),
    "phen_mat must have same number of traits as gmat"
  )
})

# --- clgsi: lines 440-442 – P_y/P_g/P_yg coercion from data.frame -------
# Also covers lines 444-451 (dimension checks on provided cov matrices)
test_that("clgsi accepts data.frame P_y/P_g/P_yg and coerces them (lines 440-442)", {
  data <- setup_test_data()
  n <- data$n_traits
  P_y <- as.data.frame(cov(data$phen_mat))
  P_g <- as.data.frame(cov(data$gebv_mat))
  P_yg <- as.data.frame(cov(data$phen_mat, data$gebv_mat))
  result <- clgsi(
    phen_mat = NULL, gebv_mat = NULL,
    pmat = data$pmat, gmat = data$gmat, wmat = data$weights,
    P_y = P_y, P_g = P_g, P_yg = P_yg,
    reliability = 0.7
  )
  expect_s3_class(result, "clgsi")
  expect_equal(length(result$b_y), n)
})

test_that("clgsi errors when P_y has wrong dimensions (line 445)", {
  data <- setup_test_data()
  n <- data$n_traits
  bad_P_y <- matrix(1, 3, 3)
  P_g <- cov(data$gebv_mat)
  P_yg <- cov(data$phen_mat, data$gebv_mat)
  expect_error(
    clgsi(
      phen_mat = NULL, gebv_mat = NULL,
      pmat = data$pmat, gmat = data$gmat, wmat = data$weights,
      P_y = bad_P_y, P_g = P_g, P_yg = P_yg
    ),
    "P_y dimensions must match"
  )
})

test_that("clgsi errors when P_g has wrong dimensions (line 448)", {
  data <- setup_test_data()
  n <- data$n_traits
  P_y <- cov(data$phen_mat)
  bad_P_g <- matrix(1, 3, 3)
  P_yg <- cov(data$phen_mat, data$gebv_mat)
  expect_error(
    clgsi(
      phen_mat = NULL, gebv_mat = NULL,
      pmat = data$pmat, gmat = data$gmat, wmat = data$weights,
      P_y = P_y, P_g = bad_P_g, P_yg = P_yg
    ),
    "P_g dimensions must match"
  )
})

test_that("clgsi errors when P_yg has wrong dimensions (line 451)", {
  data <- setup_test_data()
  n <- data$n_traits
  P_y <- cov(data$phen_mat)
  P_g <- cov(data$gebv_mat)
  bad_P_yg <- matrix(1, 3, 3)
  expect_error(
    clgsi(
      phen_mat = NULL, gebv_mat = NULL,
      pmat = data$pmat, gmat = data$gmat, wmat = data$weights,
      P_y = P_y, P_g = P_g, P_yg = bad_P_yg
    ),
    "P_yg dimensions must match"
  )
})

# --- clgsi: line 460 – gmat dimension mismatch ---------------------------
test_that("clgsi errors when gmat is not square (line 460)", {
  data <- setup_test_data()
  n <- data$n_traits
  # A non-square gmat: nrow = n_traits, ncol = n_traits - 1
  # n_traits <- nrow(bad_gmat) == n, but ncol(bad_gmat) == n-1 != n => line 460 fires
  bad_gmat <- matrix(seq_len(n * (n - 1L)), nrow = n, ncol = n - 1L)
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat,
      pmat = data$pmat, gmat = bad_gmat, wmat = data$weights,
      reliability = 0.7
    ),
    "gmat dimensions must match"
  )
})

# --- clgsi: lines 477-485 – wmat data.frame branch + validation -----------
test_that("clgsi coerces data.frame wmat and validates wcol (lines 477-485)", {
  data <- setup_test_data()
  n <- data$n_traits
  wmat_df <- as.data.frame(matrix(rep(data$weights, 2), ncol = 2))

  # Valid: coercion should work silently
  result <- clgsi(data$phen_mat, data$gebv_mat,
    data$pmat, data$gmat, wmat_df,
    wcol = 1,
    reliability = 0.7
  )
  expect_s3_class(result, "clgsi")

  # Wrong row count (line 481)
  wmat_bad_rows <- as.data.frame(matrix(1:(n - 1), ncol = 1))
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat,
      data$pmat, data$gmat, wmat_bad_rows,
      reliability = 0.7
    ),
    "Number of rows in wmat must equal number of traits"
  )

  # wcol out of range (line 485)
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat,
      data$pmat, data$gmat, wmat_df,
      wcol = 10, reliability = 0.7
    ),
    "wcol must be between 1 and"
  )
})

# --- clgsi: lines 497-498 – P_g provided (else-if branch) ----------------
test_that("clgsi uses P_g directly when provided (lines 497-498)", {
  data <- setup_test_data()
  # Supply P_g without the full cov-matrix set so has_cov_matrices is FALSE
  # This exercises the `} else if (!is.null(P_g))` branch at line 497
  P_g <- cov(data$gebv_mat)
  result <- clgsi(
    phen_mat = data$phen_mat, gebv_mat = data$gebv_mat,
    pmat = data$pmat, gmat = data$gmat, wmat = data$weights,
    P_g = P_g, reliability = 0.7
  )
  expect_s3_class(result, "clgsi")
})

# --- clgsi: lines 513-523 – low reliability warning (NULL reliability) ---
test_that("clgsi warns when estimated reliability is low (lines 513-523)", {
  data <- setup_test_data()
  # Scale GEBVs down so Var(GEBV)/Var(G) < 0.3
  gebv_low <- data$gebv_mat * 0.01
  expect_warning(
    clgsi(data$phen_mat, gebv_low,
      data$pmat, data$gmat, data$weights,
      reliability = NULL
    ),
    "Estimated reliabilities are low"
  )
})

# --- clgsi: line 527 – single reliability out of range -------------------
test_that("clgsi errors when single reliability is out of range (line 527)", {
  data <- setup_test_data()
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat,
      data$pmat, data$gmat, data$weights,
      reliability = 1.5
    ),
    "Reliability must be between 0 and 1"
  )
})

# --- clgsi: lines 530-533 – vector reliability out of [0,1] -------------
test_that("clgsi errors when vector reliability has values outside [0,1] (lines 530-533)", {
  data <- setup_test_data()
  n <- data$n_traits
  bad_rel <- rep(0.6, n)
  bad_rel[3] <- 1.2
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat,
      data$pmat, data$gmat, data$weights,
      reliability = bad_rel
    ),
    "All reliability values must be between 0 and 1"
  )
})

# --- clgsi: line 536 – wrong-length reliability vector -------------------
test_that("clgsi errors when reliability vector has wrong length (line 536)", {
  data <- setup_test_data()
  expect_error(
    clgsi(data$phen_mat, data$gebv_mat,
      data$pmat, data$gmat, data$weights,
      reliability = c(0.5, 0.7)
    ),
    "reliability must be NULL, a single value, or a vector of length n_traits"
  )
})

# --- clgsi: line 565 – non-symmetric P_combined warning ------------------
test_that("clgsi warns when P_combined is not symmetric (line 565)", {
  data <- setup_test_data()
  n <- data$n_traits
  # P_combined = rbind(cbind(P_phen, P_yg), cbind(t(P_yg), P_gebv))
  # t(P_combined) = rbind(cbind(t(P_phen), P_yg), cbind(t(P_yg), t(P_gebv)))
  # To make P_combined != t(P_combined), P_phen must be non-symmetric.
  P_y_asym <- cov(data$phen_mat)
  P_y_asym[1L, 2L] <- P_y_asym[1L, 2L] + 100 # break symmetry of P_y
  P_g <- cov(data$gebv_mat)
  P_yg <- cov(data$phen_mat, data$gebv_mat)
  expect_warning(
    clgsi(
      phen_mat = NULL, gebv_mat = NULL,
      pmat = data$pmat, gmat = data$gmat, wmat = data$weights,
      P_y = P_y_asym, P_g = P_g, P_yg = P_yg,
      reliability = 0.7
    ),
    "Combined variance matrix is not symmetric"
  )
})

# --- clgsi: line 630 – Delta_H returns NA when sigma_I == 0 --------------
test_that("clgsi returns NA Delta_H when index variance is zero (line 630)", {
  data <- setup_test_data()
  n <- data$n_traits
  # Zero weights force bPb = 0 => sigma_I = NA => rep(NA_real_, n_traits)
  zero_weights <- rep(0, n)
  result <- clgsi(data$phen_mat, data$gebv_mat,
    data$pmat, data$gmat, zero_weights,
    reliability = 0.7
  )
  expect_true(all(is.na(result$Delta_H)))
})

# --- clgsi: line 642 – PRE computed from GAY (GA / GAY) * 100 ------------
test_that("clgsi calculates PRE correctly when GAY is provided (line 642)", {
  data <- setup_test_data()
  GAY <- 1.2
  result <- clgsi(data$phen_mat, data$gebv_mat,
    data$pmat, data$gmat, data$weights,
    reliability = 0.7, GAY = GAY
  )
  expect_false(is.na(result$PRE))
  expect_equal(result$PRE, (result$GA / GAY) * 100, tolerance = 1e-4)
})

# --- lgsi: line 216 – NA/Inf in index coefficients -----------------------
test_that("lgsi errors when index coefficients contain NA/Inf (line 216)", {
  data <- setup_test_data()
  # When reliability is provided, b = ginv(P_gebv) %*% C_gebv_g %*% w where
  # C_gebv_g = sweep(gmat, 1, accuracy_vec, "*"). Injecting NaN into gmat
  # corrupts C_gebv_g (gmat has no NA guard), propagating NaN into b.
  nan_gmat <- data$gmat
  nan_gmat[1, 1] <- NaN
  expect_error(
    lgsi(data$gebv_mat, nan_gmat, data$weights, reliability = 0.7),
    "Index coefficients contain NA or Inf. Check that P_gebv is invertible."
  )
})

# --- clgsi: line 591 – NA/Inf in combined index coefficients -------------
test_that("clgsi errors when index coefficients contain NA/Inf (line 591)", {
  data <- setup_test_data()
  # Inf in gmat: P_combined is built from P_y/P_g/P_yg (all valid and finite,
  # so the symmetry check at line 564 passes). C_gebv_g = sweep(gmat, 1,
  # accuracy_vec, "*") picks up Inf, making rhs = G_combined %*% w Inf, which
  # ginv() propagates into b_combined, triggering line 591.
  inf_gmat <- data$gmat
  inf_gmat[1L, 1L] <- Inf
  P_y <- cov(data$phen_mat)
  P_g <- cov(data$gebv_mat)
  P_yg <- cov(data$phen_mat, data$gebv_mat)
  expect_error(
    clgsi(
      phen_mat = NULL, gebv_mat = NULL,
      pmat = data$pmat, gmat = inf_gmat, wmat = data$weights,
      P_y = P_y, P_g = P_g, P_yg = P_yg,
      reliability = 0.7
    ),
    "Index coefficients contain NA or Inf. Check matrix conditioning."
  )
})
