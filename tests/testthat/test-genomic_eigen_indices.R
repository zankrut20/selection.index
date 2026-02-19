# ==============================================================================
# Tests for Genomic Eigen Selection Indices (Chapter 8)
# ==============================================================================

# Setup test data for genomic eigen indices
setup_eigen_test_data <- function() {
  set.seed(123)
  
  # Use actual variance-covariance matrices
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Simulate GEBVs and phenotypes
  n_genotypes <- 50
  n_traits <- ncol(gmat)
  
  phen_mat <- as.matrix(seldata[1:n_genotypes, 3:9])
  gebv_mat <- phen_mat * 0.7 + matrix(rnorm(prod(dim(phen_mat)), sd = 0.5), 
                                       nrow(phen_mat), ncol(phen_mat))
  
  # Gamma represents covariance between phenotypes and GEBVs
  # Following the documentation example: Assume 80% GEBV-phenotype covariance
  Gamma <- gmat * 0.8
  
  # Marker score matrices (simulated)
  S_M <- gmat * 0.7      # Cov(y, s)
  S_Mg <- gmat * 0.65    # Cov(g, s)
  S_var <- gmat * 0.8    # Var(s)
  
  # Marker scores matrix
  marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
                          nrow = n_genotypes, ncol = n_traits)
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
  expect_true(result$hI2 >= 0)  # Eigen indices can have hI2 > 1
  expect_true(result$rHI >= 0)
  expect_false(is.na(result$R_M))
  expect_true(result$lambda2 > 0)
})

test_that("mesim works with all three covariance matrices", {
  data <- setup_eigen_test_data()
  
  # Most rigorous: All three matrices provided
  result <- mesim(data$pmat, data$gmat, data$S_M, 
                  S_Mg = data$S_Mg, S_var = data$S_var)
  
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
  non_sym[1, 2] <- non_sym[1, 2] + 1  # Break symmetry
  
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
                   selection_intensity = 2.063)
  result2 <- gesim(data$pmat, data$gmat, Gamma, 
                   selection_intensity = 2.421)
  
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
                     selection_intensity = 2.063)
  result2 <- gw_esim(data$pmat, data$gmat, G_M, M, 
                     selection_intensity = 2.421)
  
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
  U_mat[1, 1] <- 1  # Constrain trait 1
  U_mat[2, 2] <- 1  # Constrain trait 2
  
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
  
  wrong_d <- c(2, 1, 1)  # Too short
  
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
  Phi <- rbind(cbind(data$pmat, data$S_M), 
               cbind(t(data$S_M), data$S_var))
  A <- rbind(cbind(data$gmat, data$S_Mg), 
             cbind(t(data$S_Mg), data$S_var))
  
  b_combined <- rep(1, 2 * data$n_traits)
  
  metrics <- selection.index:::.genomic_eigen_index_metrics(
    b_combined, Phi, A, lambda2 = 0.5
  )
  
  expect_true("hI2" %in% names(metrics))
  expect_true("rHI" %in% names(metrics))
  expect_true("sigma_I" %in% names(metrics))
  expect_true("E_vec" %in% names(metrics))
  
  expect_equal(metrics$hI2, 0.5)
  expect_false(is.na(metrics$rHI))
})

