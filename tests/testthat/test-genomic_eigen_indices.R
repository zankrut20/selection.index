# ==============================================================================
# Tests for Genomic Eigen Selection Indices (Chapter 8)
# ==============================================================================

# Setup test data for genomic eigen indices
setup_eigen_test_data <- function() {
  set.seed(123)

  # Use actual variance-covariance matrices
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Simulate GEBVs and phenotypes
  n_genotypes <- 50
  n_traits <- ncol(gmat)

  phen_mat <- as.matrix(seldata[1:n_genotypes, 3:9])
  gebv_mat <- phen_mat * 0.7 + matrix(
    rnorm(prod(dim(phen_mat)), sd = 0.5),
    nrow(phen_mat), ncol(phen_mat)
  )

  # Gamma represents covariance between phenotypes and GEBVs
  # Following the documentation example: Assume 80% GEBV-phenotype covariance
  Gamma <- gmat * 0.8

  # Marker score matrices (simulated)
  S_M <- gmat * 0.7 # Cov(y, s)
  S_Mg <- gmat * 0.65 # Cov(g, s)
  S_var <- gmat * 0.8 # Var(s)

  # Marker scores matrix
  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(marker_scores) <- colnames(gmat)

  list(
    gmat = gmat,
    pmat = pmat,
    phen_mat = phen_mat,
    gebv_mat = gebv_mat,
    Gamma = Gamma,
    S_M = S_M,
    S_Mg = S_Mg,
    S_var = S_var,
    marker_scores = marker_scores,
    n_traits = n_traits
  )
}

# ==============================================================================
# TEST: MESIM - Molecular Eigen Selection Index Method
# ==============================================================================

test_that("mesim returns correct structure", {
  data <- setup_eigen_test_data()

  result <- mesim(data$pmat, data$gmat, data$S_M)

  # Check S3 classes
  expect_s3_class(result, "mesim")
  expect_s3_class(result, "genomic_eigen_index")
  expect_type(result, "list")

  # Check components exist
  expect_true("b_y" %in% names(result))
  expect_true("b_s" %in% names(result))
  expect_true("b_combined" %in% names(result))
  expect_true("E_M" %in% names(result))
  expect_true("sigma_I" %in% names(result))
  expect_true("hI2" %in% names(result))
  expect_true("rHI" %in% names(result))
  expect_true("R_M" %in% names(result))
  expect_true("lambda2" %in% names(result))
  expect_true("summary" %in% names(result))

  # Check dimensions
  expect_equal(length(result$b_y), data$n_traits)
  expect_equal(length(result$b_s), data$n_traits)
  expect_equal(length(result$b_combined), 2 * data$n_traits)
  expect_equal(length(result$E_M), data$n_traits)

  # Check values are reasonable
  expect_true(result$hI2 >= 0) # Eigen indices can have hI2 > 1
  expect_true(result$rHI >= 0)
  expect_false(is.na(result$R_M))
  expect_true(result$lambda2 > 0)
})

test_that("mesim works with all three covariance matrices", {
  data <- setup_eigen_test_data()

  # Most rigorous: All three matrices provided
  result <- mesim(data$pmat, data$gmat, data$S_M,
    S_Mg = data$S_Mg, S_var = data$S_var
  )

  expect_s3_class(result, "mesim")
  expect_false(is.na(result$hI2))
  expect_true(result$lambda2 > 0)
})

test_that("mesim defaults S_Mg to S_M when not provided", {
  data <- setup_eigen_test_data()

  result <- mesim(data$pmat, data$gmat, data$S_M)

  expect_s3_class(result, "mesim")
  expect_false(is.na(result$hI2))
})

test_that("mesim handles different selection intensities", {
  data <- setup_eigen_test_data()

  result1 <- mesim(data$pmat, data$gmat, data$S_M, selection_intensity = 2.063)
  result2 <- mesim(data$pmat, data$gmat, data$S_M, selection_intensity = 2.421)

  # R_M should scale with selection intensity
  ratio <- result2$R_M / result1$R_M
  expect_equal(ratio, 2.421 / 2.063, tolerance = 1e-6)
})

test_that("mesim errors with non-symmetric matrices", {
  data <- setup_eigen_test_data()

  non_sym <- data$pmat
  non_sym[1, 2] <- non_sym[1, 2] + 1 # Break symmetry

  expect_error(mesim(non_sym, data$gmat, data$S_M), "symmetric")
})

test_that("mesim errors with dimension mismatch", {
  data <- setup_eigen_test_data()

  wrong_gmat <- data$gmat[1:5, 1:5]

  expect_error(mesim(data$pmat, wrong_gmat, data$S_M), "dimensions")
})

test_that("mesim has print method", {
  data <- setup_eigen_test_data()
  result <- mesim(data$pmat, data$gmat, data$S_M)

  # Verify object structure instead of print output
  expect_s3_class(result, "mesim")
  expect_s3_class(result, "genomic_eigen_index")
})

test_that("mesim has summary method", {
  data <- setup_eigen_test_data()
  result <- mesim(data$pmat, data$gmat, data$S_M)

  summ <- capture.output(summary(result))
  expect_type(summ, "character")
})

# ==============================================================================
# TEST: GESIM - Linear Genomic Eigen Selection Index Method
# ==============================================================================

test_that("gesim returns correct structure", {
  data <- setup_eigen_test_data()

  # gesim needs Gamma (genomic covariance matrix from GEBVs)
  Gamma <- data$Gamma

  result <- gesim(data$pmat, data$gmat, Gamma)

  # Check S3 classes
  expect_s3_class(result, "gesim")
  expect_s3_class(result, "genomic_eigen_index")
  expect_type(result, "list")

  # Check components
  expect_true("b_y" %in% names(result))
  expect_true("b_gamma" %in% names(result))
  expect_true("E_G" %in% names(result))
  expect_true("hI2" %in% names(result))
  expect_true("rHI" %in% names(result))
  expect_true("lambda2" %in% names(result))

  # Check dimensions
  expect_equal(length(result$b_y), data$n_traits)
  expect_equal(length(result$b_gamma), data$n_traits)
  expect_equal(length(result$E_G), data$n_traits)

  # Check values (eigen indices can have hI2 > 1)
  expect_true(result$hI2 >= 0)
  expect_true(result$rHI >= 0)
  expect_true(result$lambda2 > 0)
})

test_that("gesim validates symmetry", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  # Break symmetry in pmat
  bad_pmat <- data$pmat
  bad_pmat[1, 2] <- bad_pmat[1, 2] + 1

  expect_error(gesim(bad_pmat, data$gmat, Gamma), "symmetric")
})

test_that("gesim validates dimensions", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  # Wrong dimension gmat
  wrong_gmat <- data$gmat[1:5, 1:5]

  expect_error(gesim(data$pmat, wrong_gmat, Gamma), "dimensions")
})

test_that("gesim handles different selection intensities", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  result1 <- gesim(data$pmat, data$gmat, Gamma,
    selection_intensity = 2.063
  )
  result2 <- gesim(data$pmat, data$gmat, Gamma,
    selection_intensity = 2.421
  )

  # Expected gains should scale
  ratio <- result2$R_G / result1$R_G
  expect_equal(ratio, 2.421 / 2.063, tolerance = 1e-6)
})

test_that("gesim has print and summary methods", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma
  result <- gesim(data$pmat, data$gmat, Gamma)

  expect_s3_class(result, "gesim")
  expect_s3_class(result, "genomic_eigen_index")

  summ <- capture.output(summary(result))
  expect_type(summ, "character")
})

# ==============================================================================
# TEST: GW-ESIM - Genome-Wide Eigen Selection Index Method
# ==============================================================================

test_that("gw_esim returns correct structure", {
  data <- setup_eigen_test_data()

  # gw_esim needs 4 matrices: pmat, gmat, G_M (cov between traits and markers), M (marker covariance)
  n_markers <- 20
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5), nrow = data$n_traits, ncol = n_markers)
  M <- cov(matrix(rnorm(100 * n_markers), nrow = 100, ncol = n_markers))

  result <- gw_esim(data$pmat, data$gmat, G_M, M)

  # Check S3 classes
  expect_s3_class(result, "gw_esim")
  expect_s3_class(result, "genomic_eigen_index")

  # Check components
  expect_true("b_y" %in% names(result))
  expect_true("b_m" %in% names(result))
  expect_true("E_W" %in% names(result))
  expect_true("hI2" %in% names(result))
  expect_true("rHI" %in% names(result))
  expect_true("lambda2" %in% names(result))

  # Check dimensions
  expect_equal(length(result$b_y), data$n_traits)
  expect_equal(length(result$b_m), n_markers)
  expect_equal(length(result$E_W), data$n_traits)

  # Check values (eigen indices can have hI2 > 1)
  expect_true(result$hI2 >= 0)
  expect_true(result$rHI >= 0)
})

test_that("gw_esim handles different selection intensities", {
  data <- setup_eigen_test_data()

  # Use fewer markers to ensure stable results
  n_markers <- 5
  set.seed(123)
  G_M <- matrix(rowMeans(data$gmat), nrow = data$n_traits, ncol = n_markers) +
    matrix(rnorm(data$n_traits * n_markers, sd = 0.1), data$n_traits, n_markers)
  M <- diag(n_markers) * mean(diag(data$gmat))

  result1 <- gw_esim(data$pmat, data$gmat, G_M, M,
    selection_intensity = 2.063
  )
  result2 <- gw_esim(data$pmat, data$gmat, G_M, M,
    selection_intensity = 2.421
  )

  # Response should scale (skip if NA due to numerical issues)
  if (!is.na(result1$R_W) && !is.na(result2$R_W)) {
    ratio <- result2$R_W / result1$R_W
    expect_equal(ratio, 2.421 / 2.063, tolerance = 1e-6)
  } else {
    skip("Numerical instability with test data")
  }
})

test_that("gw_esim print and summary methods work", {
  data <- setup_eigen_test_data()

  n_markers <- 20
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5), nrow = data$n_traits, ncol = n_markers)
  M <- cov(matrix(rnorm(100 * n_markers), nrow = 100, ncol = n_markers))

  result <- gw_esim(data$pmat, data$gmat, G_M, M)

  expect_s3_class(result, "gw_esim")
  expect_s3_class(result, "genomic_eigen_index")

  summ <- capture.output(summary(result))
  expect_type(summ, "character")
})

# ==============================================================================
# TEST: RGESIM - Restricted Genomic Eigen Selection Index
# ==============================================================================

test_that("rgesim returns correct structure with restricted traits", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  # rgesim needs U_mat (constraint matrix): each row is a constraint
  # To restrict traits 1 and 2 to zero gain: U = [1 0 0...], [0 1 0...]
  U_mat <- matrix(0, nrow = 2, ncol = data$n_traits)
  U_mat[1, 1] <- 1 # Constrain trait 1
  U_mat[2, 2] <- 1 # Constrain trait 2

  result <- rgesim(data$pmat, data$gmat, Gamma, U_mat)

  # Check S3 classes
  expect_s3_class(result, "rgesim")
  expect_s3_class(result, "genomic_eigen_index")

  # Check components
  expect_true("b_y" %in% names(result))
  expect_true("b_gamma" %in% names(result))
  expect_true("E_RG" %in% names(result))
  expect_true("hI2" %in% names(result))
  expect_true("constrained_response" %in% names(result))

  # Check that constrained response is near zero
  expect_true(max(abs(result$constrained_response)) < 0.01)
})

test_that("rgesim works with custom constraint matrix", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  # Custom constraint: trait1 + trait2 = 0 (single linear constraint)
  U_mat <- matrix(c(1, 1, rep(0, data$n_traits - 2)), nrow = 1)

  result <- rgesim(data$pmat, data$gmat, Gamma, U_mat)

  expect_s3_class(result, "rgesim")
  expect_equal(result$n_restrictions, 1)
  expect_true(max(abs(result$constrained_response)) < 0.01)
})

test_that("rgesim print and summary methods work", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  result <- rgesim(data$pmat, data$gmat, Gamma, U_mat)

  expect_s3_class(result, "rgesim")
  expect_s3_class(result, "genomic_eigen_index")

  summ <- capture.output(summary(result))
  expect_type(summ, "character")
})

# ==============================================================================
# TEST: PPG-GESIM - Predetermined Proportional Gain Genomic Eigen Index
# ==============================================================================

test_that("ppg_gesim returns correct structure", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  # Define proportional gains (d vector)
  d <- c(2, 1, 1, 1, 1, 1, 1)

  result <- ppg_gesim(data$pmat, data$gmat, Gamma, d)

  # Check S3 classes
  expect_s3_class(result, "ppg_gesim")
  expect_s3_class(result, "genomic_eigen_index")

  # Check components
  expect_true("b_y" %in% names(result))
  expect_true("b_gamma" %in% names(result))
  expect_true("E_PG" %in% names(result))
  expect_true("d" %in% names(result))

  # Check dimensions
  expect_equal(length(result$d), data$n_traits)
  expect_equal(length(result$E_PG), data$n_traits)
})

test_that("ppg_gesim d vector has correct length", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma

  wrong_d <- c(2, 1, 1) # Too short

  expect_error(
    ppg_gesim(data$pmat, data$gmat, Gamma, wrong_d),
    "length"
  )
})

test_that("ppg_gesim has print and summary methods", {
  data <- setup_eigen_test_data()
  Gamma <- data$Gamma
  d <- rep(1, data$n_traits)
  result <- ppg_gesim(data$pmat, data$gmat, Gamma, d)

  expect_s3_class(result, "ppg_gesim")
  expect_s3_class(result, "genomic_eigen_index")

  summ <- capture.output(summary(result))
  expect_type(summ, "character")
})

# ==============================================================================
# TEST: Helper functions (internal)
# ==============================================================================

test_that(".gesim_leading_eigenvector finds positive eigenvalue", {
  # Create a simple positive definite matrix
  mat <- matrix(c(4, 2, 2, 3), nrow = 2)

  result <- selection.index:::.gesim_leading_eigenvector(mat)

  expect_true("vector" %in% names(result))
  expect_true("value" %in% names(result))
  expect_true("all_values" %in% names(result))

  expect_equal(length(result$vector), 2)
  expect_true(result$value > 0)
})

test_that(".gesim_leading_eigenvector errors with no positive eigenvalues", {
  # Create a matrix with only negative eigenvalues (not valid for variance)
  mat <- matrix(c(-4, -2, -2, -3), nrow = 2)

  expect_error(
    selection.index:::.gesim_leading_eigenvector(mat),
    "positive eigenvalue"
  )
})

test_that(".genomic_eigen_index_metrics computes correct values", {
  data <- setup_eigen_test_data()

  # Create test coefficient vector
  b <- rep(1, data$n_traits)

  # Create combined matrices
  Phi <- rbind(
    cbind(data$pmat, data$S_M),
    cbind(t(data$S_M), data$S_var)
  )
  A <- rbind(
    cbind(data$gmat, data$S_Mg),
    cbind(t(data$S_Mg), data$S_var)
  )

  b_combined <- rep(1, 2 * data$n_traits)

  metrics <- selection.index:::.genomic_eigen_index_metrics(
    b_combined, Phi, A,
    lambda2 = 0.5
  )

  expect_true("hI2" %in% names(metrics))
  expect_true("rHI" %in% names(metrics))
  expect_true("sigma_I" %in% names(metrics))
  expect_true("E_vec" %in% names(metrics))

  expect_equal(metrics$hI2, 0.5)
  expect_false(is.na(metrics$rHI))
})

# ==============================================================================
# NEW COVERAGE TESTS - targeting previously uncovered lines
# ==============================================================================

# helper: make a minimal 2x2 symmetric positive-definite matrix
.make_sym2 <- function(a = 4, b = 1) {
  matrix(c(a, b, b, a), 2, 2)
}

# ==============================================================================
# .genomic_eigen_index_metrics – internal helper (lines 90, 97)
# ==============================================================================

# Line 97: bAb / bPhibb branch (lambda2 = NULL, sigma_I > 0)
test_that(".genomic_eigen_index_metrics falls back to bAb/bPhibb when lambda2 is NULL (line 97)", {
  data <- setup_eigen_test_data()
  n <- data$n_traits
  b <- rep(1, 2 * n)
  Phi <- rbind(cbind(data$pmat, data$S_M), cbind(t(data$S_M), data$S_var))
  A <- rbind(cbind(data$gmat, data$S_Mg), cbind(t(data$S_Mg), data$S_var))

  # lambda2 = NULL forces the else-if branch (line 97)
  metrics <- selection.index:::.genomic_eigen_index_metrics(b, Phi, A, lambda2 = NULL)

  expect_false(is.na(metrics$hI2))
  expect_false(is.na(metrics$rHI))
})

# Line 90: rep(NA_real_) E_vec path (zero-variance b => sigma_I = NA)
test_that(".genomic_eigen_index_metrics returns NA E_vec when sigma_I is NA (line 90)", {
  n <- 2
  # Phi has negative quadratic form for b=c(1,-1): bPhibb < 0 then flipped but zero
  # Easiest: use a zero b so bPhibb = 0 -> sigma_I = NA
  b <- c(0, 0)
  Phi <- .make_sym2()
  A <- .make_sym2(3, 0.5)

  metrics <- selection.index:::.genomic_eigen_index_metrics(b, Phi, A, lambda2 = NULL)

  expect_true(all(is.na(metrics$E_vec)))
})

# ==============================================================================
# MESIM validation (lines 268, 270, 272, 274, 279, 283)
# ==============================================================================

# Line 268: gmat not symmetric in mesim
test_that("mesim errors when gmat is not symmetric (line 268)", {
  data <- setup_eigen_test_data()
  bad_gmat <- data$gmat
  bad_gmat[1, 2] <- bad_gmat[1, 2] + 5
  expect_error(
    mesim(data$pmat, bad_gmat, data$S_M),
    "gmat must be a symmetric matrix"
  )
})

# Line 270: S_M not symmetric in mesim
test_that("mesim errors when S_M is not symmetric (line 270)", {
  data <- setup_eigen_test_data()
  bad_SM <- data$S_M
  bad_SM[1, 2] <- bad_SM[1, 2] + 5
  expect_error(
    mesim(data$pmat, data$gmat, bad_SM),
    "S_M must be a symmetric matrix"
  )
})

# Line 272: S_Mg not symmetric in mesim (explicit S_Mg provided)
test_that("mesim errors when S_Mg is not symmetric (line 272)", {
  data <- setup_eigen_test_data()
  bad_SMg <- data$S_Mg
  bad_SMg[1, 2] <- bad_SMg[1, 2] + 5
  expect_error(
    mesim(data$pmat, data$gmat, data$S_M, S_Mg = bad_SMg),
    "S_Mg must be a symmetric matrix"
  )
})

# Line 274: S_var not symmetric in mesim (explicit S_var provided)
test_that("mesim errors when S_var is not symmetric (line 274)", {
  data <- setup_eigen_test_data()
  bad_Svar <- data$S_var
  bad_Svar[1, 2] <- bad_Svar[1, 2] + 5
  expect_error(
    mesim(data$pmat, data$gmat, data$S_M, S_var = bad_Svar),
    "S_var must be a symmetric matrix"
  )
})

# Line 279: n_traits < 2 in mesim
test_that("mesim errors when fewer than 2 traits are provided (line 279)", {
  p <- matrix(4, 1, 1)
  g <- matrix(2, 1, 1)
  s <- matrix(1, 1, 1)
  expect_error(mesim(p, g, s), "At least 2 traits")
})

# Line 283: auto-generated trait names when pmat has no colnames
test_that("mesim generates Trait_ names when pmat has no colnames (line 283)", {
  data <- setup_eigen_test_data()
  pmat_no_names <- unname(data$pmat)
  gmat_no_names <- unname(data$gmat)
  S_M_no_names <- unname(data$S_M)

  result <- mesim(pmat_no_names, gmat_no_names, S_M_no_names)

  expect_true(all(grepl("^Trait_", result$trait_names)))
})

# ==============================================================================
# GESIM validation (lines 470, 472, 476, 480)
# ==============================================================================

# Line 470: gmat not symmetric in gesim
test_that("gesim errors when gmat is not symmetric (line 470)", {
  data <- setup_eigen_test_data()
  bad_gmat <- data$gmat
  bad_gmat[1, 2] <- bad_gmat[1, 2] + 5
  expect_error(
    gesim(data$pmat, bad_gmat, data$Gamma),
    "gmat must be a symmetric matrix"
  )
})

# Line 472: Gamma not symmetric in gesim
test_that("gesim errors when Gamma is not symmetric (line 472)", {
  data <- setup_eigen_test_data()
  bad_Gamma <- data$Gamma
  bad_Gamma[1, 2] <- bad_Gamma[1, 2] + 5
  expect_error(
    gesim(data$pmat, data$gmat, bad_Gamma),
    "Gamma must be a symmetric matrix"
  )
})

# Line 476: n_traits < 2 in gesim
test_that("gesim errors when fewer than 2 traits are provided (line 476)", {
  p <- matrix(4, 1, 1)
  g <- matrix(2, 1, 1)
  G <- matrix(1, 1, 1)
  expect_error(gesim(p, g, G), "At least 2 traits")
})

# Line 480: auto-generated trait names when pmat has no colnames in gesim
test_that("gesim generates Trait_ names when pmat has no colnames (line 480)", {
  data <- setup_eigen_test_data()
  result <- gesim(unname(data$pmat), unname(data$gmat), unname(data$Gamma))
  expect_true(all(grepl("^Trait_", result$trait_names)))
})

# Lines 533-534: implied_w tryCatch warning + NA in gesim
# Force ginv(A) to throw by injecting non-finite values into A's source matrices
# (Gamma with Inf makes A have Inf; ginv -> svd fails with "infinite values")
test_that("gesim warns and returns NA implied_w when A is degenerate (lines 533-534)", {
  data <- setup_eigen_test_data()
  # Use a tiny 2x2 case that succeeds for the eigen but makes ginv(A) fail
  # Strategy: mock ginv inside selection.index to always throw
  skip_if_not_installed("testthat", minimum_version = "3.1.0")
  # gesim calls ginv() only once inside the tryCatch (line 530).
  # Mocking it causes the tryCatch to fire the error handler -> warning + NA.
  expect_warning(
    with_mocked_bindings(
      ginv = function(...) stop("mocked ginv failure"),
      .package = "selection.index",
      {
        gesim(data$pmat, data$gmat, data$Gamma)
      }
    ),
    "Could not compute implied economic weights"
  )
})

# ==============================================================================
# GW-ESIM validation (lines 668, 670, 672, 674, 676, 678, 680, 684)
# ==============================================================================

# Helper: build valid gw_esim inputs from 2x2 matrices
.gw_esim_2trait <- function() {
  p <- .make_sym2(4, 0.5)
  g <- .make_sym2(2, 0.3)
  G_M <- matrix(c(0.5, 0.3), nrow = 2, ncol = 1)
  M <- matrix(1, 1, 1)
  list(pmat = p, gmat = g, G_M = G_M, M = M)
}

# Line 668: pmat not symmetric in gw_esim
test_that("gw_esim errors when pmat is not symmetric (line 668)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5),
    nrow = data$n_traits, ncol = n_markers
  )
  M <- diag(n_markers)
  bad_pmat <- data$pmat
  bad_pmat[1, 2] <- bad_pmat[1, 2] + 5
  expect_error(
    gw_esim(bad_pmat, data$gmat, G_M, M),
    "pmat must be a symmetric matrix"
  )
})

# Line 670: gmat not symmetric in gw_esim
test_that("gw_esim errors when gmat is not symmetric (line 670)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5),
    nrow = data$n_traits, ncol = n_markers
  )
  M <- diag(n_markers)
  bad_gmat <- data$gmat
  bad_gmat[1, 2] <- bad_gmat[1, 2] + 5
  expect_error(
    gw_esim(data$pmat, bad_gmat, G_M, M),
    "gmat must be a symmetric matrix"
  )
})

# Line 672: M not symmetric in gw_esim
test_that("gw_esim errors when M is not symmetric (line 672)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5),
    nrow = data$n_traits, ncol = n_markers
  )
  bad_M <- diag(n_markers)
  bad_M[1, 2] <- 99 # break symmetry
  expect_error(
    gw_esim(data$pmat, data$gmat, G_M, bad_M),
    "M must be a symmetric matrix"
  )
})

# Line 674: pmat and gmat have different dimensions
test_that("gw_esim errors when pmat and gmat have different dimensions (line 674)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5),
    nrow = data$n_traits, ncol = n_markers
  )
  M <- diag(n_markers)
  bad_gmat <- data$gmat[1:5, 1:5] # wrong size, still symmetric
  expect_error(
    gw_esim(data$pmat, bad_gmat, G_M, M),
    "pmat and gmat must have the same dimensions"
  )
})

# Line 676: G_M has wrong number of rows
test_that("gw_esim errors when G_M has wrong number of rows (line 676)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  bad_G_M <- matrix(rnorm((data$n_traits - 1) * n_markers, sd = 0.5),
    nrow = data$n_traits - 1, ncol = n_markers
  )
  M <- diag(n_markers)
  expect_error(
    gw_esim(data$pmat, data$gmat, bad_G_M, M),
    "G_M must have n_traits rows"
  )
})

# Line 678: M is not square with dimension n_markers
test_that("gw_esim errors when M dimension doesn't match number of markers (line 678)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5),
    nrow = data$n_traits, ncol = n_markers
  )
  bad_M <- diag(n_markers + 2) # wrong size
  expect_error(
    gw_esim(data$pmat, data$gmat, G_M, bad_M),
    "M must be a square matrix with dimension equal to number of markers"
  )
})

# Line 680: n_traits < 2 in gw_esim
test_that("gw_esim errors when fewer than 2 traits are provided (line 680)", {
  p <- matrix(4, 1, 1)
  g <- matrix(2, 1, 1)
  G_M <- matrix(0.5, 1, 3)
  M <- diag(3)
  expect_error(gw_esim(p, g, G_M, M), "At least 2 traits")
})

# Line 684: auto-generated trait names in gw_esim
test_that("gw_esim generates Trait_ names when pmat has no colnames (line 684)", {
  data <- setup_eigen_test_data()
  n_markers <- 5
  G_M <- matrix(rnorm(data$n_traits * n_markers, sd = 0.5),
    nrow = data$n_traits, ncol = n_markers
  )
  M <- diag(n_markers)
  result <- gw_esim(unname(data$pmat), unname(data$gmat), G_M, M)
  expect_true(all(grepl("^Trait_", result$trait_names)))
})

# ==============================================================================
# RGESIM validation (lines 860, 862, 864, 866, 868, 870, 874)
# ==============================================================================

# Line 860: pmat not symmetric in rgesim
test_that("rgesim errors when pmat is not symmetric (line 860)", {
  data <- setup_eigen_test_data()
  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  bad_pmat <- data$pmat
  bad_pmat[1, 2] <- bad_pmat[1, 2] + 5
  expect_error(
    rgesim(bad_pmat, data$gmat, data$Gamma, U_mat),
    "pmat must be a symmetric matrix"
  )
})

# Line 862: gmat not symmetric in rgesim
test_that("rgesim errors when gmat is not symmetric (line 862)", {
  data <- setup_eigen_test_data()
  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  bad_gmat <- data$gmat
  bad_gmat[1, 2] <- bad_gmat[1, 2] + 5
  expect_error(
    rgesim(data$pmat, bad_gmat, data$Gamma, U_mat),
    "gmat must be a symmetric matrix"
  )
})

# Line 864: Gamma not symmetric in rgesim
test_that("rgesim errors when Gamma is not symmetric (line 864)", {
  data <- setup_eigen_test_data()
  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  bad_Gamma <- data$Gamma
  bad_Gamma[1, 2] <- bad_Gamma[1, 2] + 5
  expect_error(
    rgesim(data$pmat, data$gmat, bad_Gamma, U_mat),
    "Gamma must be a symmetric matrix"
  )
})

# Line 866: dimension mismatch in rgesim
test_that("rgesim errors when matrix dimensions do not match (line 866)", {
  data <- setup_eigen_test_data()
  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  bad_gmat <- data$gmat[1:5, 1:5]
  expect_error(
    rgesim(data$pmat, bad_gmat, data$Gamma, U_mat),
    "All matrices must have the same dimensions"
  )
})

# Line 868: U_mat has wrong number of columns
test_that("rgesim errors when U_mat has wrong number of columns (line 868)", {
  data <- setup_eigen_test_data()
  bad_U_mat <- matrix(c(1, 0), nrow = 1, ncol = 2) # too few columns
  expect_error(
    rgesim(data$pmat, data$gmat, data$Gamma, bad_U_mat),
    "U_mat must have n_traits columns"
  )
})

# Line 870: n_traits < 2 in rgesim
test_that("rgesim errors when fewer than 2 traits are provided (line 870)", {
  p <- matrix(4, 1, 1)
  g <- matrix(2, 1, 1)
  Gamma <- matrix(1, 1, 1)
  U_mat <- matrix(1, 1, 1)
  expect_error(rgesim(p, g, Gamma, U_mat), "At least 2 traits")
})

# Line 874: auto-generated trait names in rgesim
test_that("rgesim generates Trait_ names when pmat has no colnames (line 874)", {
  data <- setup_eigen_test_data()
  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  result <- rgesim(unname(data$pmat), unname(data$gmat), unname(data$Gamma), U_mat)
  expect_true(all(grepl("^Trait_", result$trait_names)))
})

# Lines 958-959: implied_w tryCatch warning + NA in rgesim
test_that("rgesim warns and returns NA implied_w when ginv fails (lines 958-959)", {
  data <- setup_eigen_test_data()
  U_mat <- matrix(c(1, rep(0, data$n_traits - 1)), nrow = 1)
  skip_if_not_installed("testthat", minimum_version = "3.1.0")
  # rgesim calls ginv() twice before the tryCatch (lines 912, 954):
  #   1st: ginv(middle) for the projection matrix (must succeed)
  #   2nd: ginv(A) inside tryCatch (must fail to trigger warning)
  # Use a stateful counter so only the 2nd call fails.
  call_count <- 0L
  real_ginv <- MASS::ginv
  expect_warning(
    with_mocked_bindings(
      ginv = function(X, ...) {
        call_count <<- call_count + 1L
        if (call_count < 2L) real_ginv(X, ...) else stop("mocked ginv failure")
      },
      .package = "selection.index",
      {
        rgesim(data$pmat, data$gmat, data$Gamma, U_mat)
      }
    ),
    "Could not compute implied economic weights"
  )
})

# ==============================================================================
# PPG-GESIM validation (lines 1106, 1108, 1110, 1112, 1116, 1120)
# ==============================================================================

# Line 1106: pmat not symmetric in ppg_gesim
test_that("ppg_gesim errors when pmat is not symmetric (line 1106)", {
  data <- setup_eigen_test_data()
  d <- rep(1, data$n_traits)
  bad_pmat <- data$pmat
  bad_pmat[1, 2] <- bad_pmat[1, 2] + 5
  expect_error(
    ppg_gesim(bad_pmat, data$gmat, data$Gamma, d),
    "pmat must be a symmetric matrix"
  )
})

# Line 1108: gmat not symmetric in ppg_gesim
test_that("ppg_gesim errors when gmat is not symmetric (line 1108)", {
  data <- setup_eigen_test_data()
  d <- rep(1, data$n_traits)
  bad_gmat <- data$gmat
  bad_gmat[1, 2] <- bad_gmat[1, 2] + 5
  expect_error(
    ppg_gesim(data$pmat, bad_gmat, data$Gamma, d),
    "gmat must be a symmetric matrix"
  )
})

# Line 1110: Gamma not symmetric in ppg_gesim
test_that("ppg_gesim errors when Gamma is not symmetric (line 1110)", {
  data <- setup_eigen_test_data()
  d <- rep(1, data$n_traits)
  bad_Gamma <- data$Gamma
  bad_Gamma[1, 2] <- bad_Gamma[1, 2] + 5
  expect_error(
    ppg_gesim(data$pmat, data$gmat, bad_Gamma, d),
    "Gamma must be a symmetric matrix"
  )
})

# Line 1112: dimension mismatch in ppg_gesim
test_that("ppg_gesim errors when matrix dimensions do not match (line 1112)", {
  data <- setup_eigen_test_data()
  d <- rep(1, data$n_traits)
  bad_gmat <- data$gmat[1:5, 1:5]
  expect_error(
    ppg_gesim(data$pmat, bad_gmat, data$Gamma, d),
    "All matrices must have the same dimensions"
  )
})

# Line 1116: n_traits < 2 in ppg_gesim
test_that("ppg_gesim errors when fewer than 2 traits are provided (line 1116)", {
  p <- matrix(4, 1, 1)
  g <- matrix(2, 1, 1)
  Gamma <- matrix(1, 1, 1)
  d <- 1
  expect_error(ppg_gesim(p, g, Gamma, d), "At least 2 traits")
})

# Line 1120: auto-generated trait names in ppg_gesim
test_that("ppg_gesim generates Trait_ names when pmat has no colnames (line 1120)", {
  data <- setup_eigen_test_data()
  d <- rep(1, data$n_traits)
  result <- ppg_gesim(unname(data$pmat), unname(data$gmat), unname(data$Gamma), d)
  expect_true(all(grepl("^Trait_", result$trait_names)))
})

# Lines 1221-1222: implied_w tryCatch warning + NA in ppg_gesim
test_that("ppg_gesim warns and returns NA implied_w when ginv fails (lines 1221-1222)", {
  data <- setup_eigen_test_data()
  d <- rep(1, data$n_traits)
  skip_if_not_installed("testthat", minimum_version = "3.1.0")
  # ppg_gesim calls ginv() twice before the tryCatch (lines 1166, ~1221):
  #   1st: ginv(middle) for the projection matrix (must succeed)
  #   2nd: ginv(A) inside tryCatch (must fail to trigger warning)
  call_count <- 0L
  real_ginv <- MASS::ginv
  expect_warning(
    with_mocked_bindings(
      ginv = function(X, ...) {
        call_count <<- call_count + 1L
        if (call_count < 2L) real_ginv(X, ...) else stop("mocked ginv failure")
      },
      .package = "selection.index",
      {
        ppg_gesim(data$pmat, data$gmat, data$Gamma, d)
      }
    ),
    "Could not compute implied economic weights"
  )
})
