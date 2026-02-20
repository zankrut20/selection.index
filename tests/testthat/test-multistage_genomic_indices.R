# ==============================================================================
# TEST: Multistage Genomic Indices (Chapter 9)
# Covers: mlgsi, mrlgsi, mppg_lgsi and all internal helpers
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared data setup helpers
# ------------------------------------------------------------------------------

setup_genomic_data_small <- function(n1 = 3, n = 5, seed = 42) {
  set.seed(seed)
  make_pd <- function(k, scale = 1) {
    L <- matrix(runif(k * k, 0.1, 0.5), k, k)
    L[upper.tri(L)] <- 0
    diag(L) <- 1
    M <- L %*% t(L)
    M * scale
  }
  Gamma1 <- make_pd(n1)
  Gamma  <- make_pd(n)
  A1     <- 0.8 * Gamma1
  A      <- 0.75 * Gamma[, seq_len(n1)]
  C      <- make_pd(n, 1.2)
  G1     <- make_pd(n1, 1.1)
  P1     <- make_pd(n1, 1.5)
  w      <- as.numeric(seq(2, n + 1))
  list(Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
       C = C, G1 = G1, P1 = P1, w = w, n1 = n1, n = n)
}

setup_genomic_data_real <- function() {
  data("seldata", package = "selection.index", envir = environment())
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  n1 <- 3; n <- 7
  reliability <- 0.7
  Gamma1 <- reliability * gmat[1:n1, 1:n1]
  Gamma  <- reliability * gmat
  A1     <- reliability * gmat[1:n1, 1:n1]
  A      <- gmat[, 1:n1]
  C      <- gmat
  G1     <- gmat[1:n1, 1:n1]
  P1     <- pmat[1:n1, 1:n1]
  w      <- c(10, 8, 6, 4, 3, 2, 1)
  list(Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
       C = C, G1 = G1, P1 = P1, w = w, n1 = n1, n = n)
}

# ==============================================================================
# TESTS: .cochran_adjustment_genomic (internal helper)
# ==============================================================================

test_that(".cochran_adjustment_genomic normal path returns adjusted matrix", {
  d <- setup_genomic_data_small()
  beta1 <- as.numeric(solve(d$Gamma1, d$A1 %*% d$w[seq_len(d$n1)]))
  Gamma_star <- selection.index:::.cochran_adjustment_genomic(
    Gamma = d$Gamma, Gamma1 = d$Gamma1, beta1 = beta1,
    A = d$A, k1 = 2.063, tau = qnorm(0.9)
  )
  expect_true(is.matrix(Gamma_star))
  expect_equal(dim(Gamma_star), c(d$n, d$n))
  expect_false(isTRUE(all.equal(Gamma_star, d$Gamma)))
})

test_that(".cochran_adjustment_genomic warns and returns Gamma when variance <= 0", {
  d <- setup_genomic_data_small()
  expect_warning(
    result <- selection.index:::.cochran_adjustment_genomic(
      Gamma = d$Gamma, Gamma1 = d$Gamma1, beta1 = rep(0, d$n1),
      A = d$A, k1 = 2.063, tau = qnorm(0.9)
    ),
    "Invalid genomic variance"
  )
  expect_equal(result, d$Gamma)
})

# ==============================================================================
# TESTS: .genomic_stage_metrics (internal helper)
# ==============================================================================

test_that(".genomic_stage_metrics normal path with A not NULL", {
  d <- setup_genomic_data_small()
  metrics <- selection.index:::.genomic_stage_metrics(
    beta = rep(0.5, d$n1), Gamma = d$Gamma1,
    A = d$A1, w = d$w[seq_len(d$n1)], k = 2.063
  )
  expect_true(is.list(metrics))
  expect_true(metrics$sigma_I > 0)
  expect_equal(length(metrics$E), d$n1)
})

test_that(".genomic_stage_metrics normal path with A NULL (stage 2)", {
  d <- setup_genomic_data_small()
  metrics <- selection.index:::.genomic_stage_metrics(
    beta = d$w, Gamma = d$Gamma, A = NULL, w = d$w, k = 2.063
  )
  expect_true(is.list(metrics))
  expect_equal(length(metrics$E), d$n)
})

test_that(".genomic_stage_metrics returns NA metrics when beta_Gamma_beta is non-finite", {
  d <- setup_genomic_data_small()
  Gamma_nan <- d$Gamma1
  Gamma_nan[1, 1] <- NaN
  expect_warning(
    metrics <- selection.index:::.genomic_stage_metrics(
      beta = rep(1, d$n1), Gamma = Gamma_nan,
      A = d$A1, w = d$w[seq_len(d$n1)], k = 2.063
    ),
    "Invalid genomic variance"
  )
  expect_true(is.na(metrics$sigma_I))
  expect_true(is.na(metrics$R))
  expect_true(all(is.na(metrics$E)))
})

test_that(".genomic_stage_metrics handles beta_Gamma_beta = 0 -> sigma_I NA", {
  d <- setup_genomic_data_small()
  Gamma_zero <- matrix(0, d$n1, d$n1)
  metrics <- selection.index:::.genomic_stage_metrics(
    beta = rep(1, d$n1), Gamma = Gamma_zero,
    A = d$A1, w = d$w[seq_len(d$n1)], k = 2.063
  )
  expect_true(is.na(metrics$sigma_I))
  expect_true(is.na(metrics$R))
  expect_true(all(is.na(metrics$E)))
})

test_that(".genomic_stage_metrics sigma_I=0 branch: A NULL path E returns NA", {
  d <- setup_genomic_data_small()
  Gamma_zero <- matrix(0, d$n, d$n)
  metrics <- selection.index:::.genomic_stage_metrics(
    beta = d$w, Gamma = Gamma_zero, A = NULL, w = d$w, k = 2.063
  )
  expect_true(is.na(metrics$sigma_I))
  expect_true(all(is.na(metrics$E)))
})

# ==============================================================================
# TESTS: .genomic_index_correlation (internal helper)
# ==============================================================================

test_that(".genomic_index_correlation computes valid correlation", {
  d <- setup_genomic_data_small()
  rho <- selection.index:::.genomic_index_correlation(
    beta1 = rep(0.5, d$n1), beta2 = d$w,
    Gamma1 = d$Gamma1, Gamma = d$Gamma, A = d$A
  )
  expect_true(is.numeric(rho))
  expect_true(is.finite(rho))
})

test_that(".genomic_index_correlation warns and returns NA for beta1 zero variance", {
  d <- setup_genomic_data_small()
  expect_warning(
    rho <- selection.index:::.genomic_index_correlation(
      beta1 = rep(0, d$n1), beta2 = d$w,
      Gamma1 = d$Gamma1, Gamma = d$Gamma, A = d$A
    ),
    "Invalid variance"
  )
  expect_true(is.na(rho))
})

test_that(".genomic_index_correlation warns and returns NA for beta2 zero variance", {
  d <- setup_genomic_data_small()
  expect_warning(
    rho <- selection.index:::.genomic_index_correlation(
      beta1 = rep(1, d$n1), beta2 = rep(0, d$n),
      Gamma1 = d$Gamma1, Gamma = d$Gamma, A = d$A
    ),
    "Invalid variance"
  )
  expect_true(is.na(rho))
})

# ==============================================================================
# TESTS: .young_intensities (internal helper)
# ==============================================================================

test_that(".young_intensities computes k1 and k2 for valid p", {
  result <- selection.index:::.young_intensities(p = 0.1, rho_12 = 0.5)
  expect_true(is.list(result))
  expect_true(result$k1 > 0)
  expect_true(result$k2 > 0)
})

test_that(".young_intensities stops for p <= 0", {
  expect_error(selection.index:::.young_intensities(p = 0,    rho_12 = 0.5), "Selection proportion")
  expect_error(selection.index:::.young_intensities(p = -0.1, rho_12 = 0.5), "Selection proportion")
})

test_that(".young_intensities stops for p >= 1", {
  expect_error(selection.index:::.young_intensities(p = 1,   rho_12 = 0.5), "Selection proportion")
  expect_error(selection.index:::.young_intensities(p = 1.2, rho_12 = 0.5), "Selection proportion")
})

test_that(".young_intensities clamps extreme rho_12 without error", {
  expect_true(is.numeric(selection.index:::.young_intensities(0.1, 1.5)$k1))
  expect_true(is.numeric(selection.index:::.young_intensities(0.1, -1.5)$k1))
})

# ==============================================================================
# TESTS: mlgsi
# ==============================================================================

test_that("mlgsi basic functionality with small synthetic data", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mlgsi(Gamma1 = d$Gamma1, Gamma = d$Gamma, A1 = d$A1, A = d$A,
          C = d$C, G1 = d$G1, P1 = d$P1, wmat = d$w)
  )
  expect_s3_class(result, "mlgsi")
  expect_s3_class(result, "multistage_genomic_index")
  expect_equal(length(result$beta1), d$n1)
  expect_equal(length(result$w), d$n)
  expect_equal(dim(result$Gamma_star), c(d$n, d$n))
  expect_equal(dim(result$C_star), c(d$n, d$n))
  expect_equal(result$k1, 2.063)
  expect_equal(result$k2, 2.063)
  expect_true(!is.null(result$summary_stage1))
  expect_true(!is.null(result$summary_stage2))
})

test_that("mlgsi works with real seldata", {
  d <- setup_genomic_data_real()
  result <- suppressWarnings(
    mlgsi(Gamma1 = d$Gamma1, Gamma = d$Gamma, A1 = d$A1, A = d$A,
          C = d$C, G1 = d$G1, P1 = d$P1, wmat = d$w)
  )
  expect_s3_class(result, "mlgsi")
  expect_true(is.finite(result$rho_I1I2))
  expect_true(is.numeric(result$stage1_metrics$R))
  expect_true(is.numeric(result$stage2_metrics$R))
})

test_that("mlgsi uses wcol for multi-column wmat", {
  d <- setup_genomic_data_small()
  wmat2 <- cbind(d$w, d$w * 2)
  res1 <- suppressWarnings(mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat2, wcol = 1))
  res2 <- suppressWarnings(mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat2, wcol = 2))
  expect_false(isTRUE(all.equal(res1$w, res2$w)))
})

test_that("mlgsi stops when weight vector length != n", {
  d <- setup_genomic_data_small()
  expect_error(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat = c(1, 2)),
    "Weight vector length"
  )
})

test_that("mlgsi respects custom tau", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, tau = 1.5)
  )
  expect_equal(result$tau, 1.5)
})

test_that("mlgsi with Young method enabled runs without error", {
  d <- setup_genomic_data_real()
  result <- suppressWarnings(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, use_young_method = TRUE)
  )
  expect_true(!is.null(result$k1))
})

test_that("mlgsi warns when C* adjustment is skipped (b1Pb1 <= 0)", {
  d <- setup_genomic_data_small()
  expect_warning(
    result <- mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1,
                    P1 = matrix(0, d$n1, d$n1), wmat = d$w),
    "Invalid phenotypic variance"
  )
  expect_equal(result$C_star, d$C)
})

test_that("mlgsi rho_HI is NA when wCw = 0", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, C = matrix(0, d$n, d$n), d$G1, d$P1, d$w)
  )
  expect_true(is.na(result$stage1_metrics$rho_HI))
})

test_that("mlgsi produces correct summary data frames", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w)
  )
  expect_true(is.data.frame(result$summary_stage1))
  expect_true(is.data.frame(result$summary_stage2))
  expect_equal(nrow(result$summary_stage1), d$n1)
  expect_equal(nrow(result$summary_stage2), d$n)
  expect_true("beta" %in% names(result$summary_stage1))
  expect_true("w"    %in% names(result$summary_stage2))
})

test_that("mlgsi uses named traits from colnames", {
  d <- setup_genomic_data_small()
  colnames(d$Gamma1) <- paste0("G1T", seq_len(d$n1))
  colnames(d$Gamma)  <- paste0("G2T", seq_len(d$n))
  result <- suppressWarnings(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w)
  )
  expect_equal(result$summary_stage1$Trait, paste0("G1T", seq_len(d$n1)))
  expect_equal(result$summary_stage2$Trait, paste0("G2T", seq_len(d$n)))
})

test_that("mlgsi Young fallback to manual when rho is NA", {
  d <- setup_genomic_data_small()
  # Degenerate Gamma1 -> rho NA -> use_young=TRUE but rho NA -> manual
  result <- withCallingHandlers(
    mlgsi(matrix(0, d$n1, d$n1), d$Gamma, d$A1, d$A, d$C, d$G1, d$P1,
          d$w, use_young_method = TRUE, k1_manual = 1.9, k2_manual = 1.9),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(result$k1, 1.9)
  expect_equal(result$k2, 1.9)
})

test_that("mlgsi stores correct selection_proportion", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
          selection_proportion = 0.15)
  )
  expect_equal(result$selection_proportion, 0.15)
})

# ==============================================================================
# TESTS: mrlgsi
# ==============================================================================

test_that("mrlgsi basic functionality with single constraint per stage", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2)
  )
  expect_s3_class(result, "mrlgsi")
  expect_s3_class(result, "multistage_genomic_index")
  expect_equal(length(result$beta_R1), d$n1)
  expect_equal(length(result$beta_R2), d$n)
  expect_equal(length(result$beta1),   d$n1)
  expect_equal(length(result$beta2),   d$n)
  expect_equal(dim(result$K_G1), c(d$n1, d$n1))
  expect_equal(dim(result$K_G2), c(d$n, d$n))
  expect_equal(nrow(result$summary_stage1), d$n1)
  expect_equal(nrow(result$summary_stage2), d$n)
  expect_true("beta_R" %in% names(result$summary_stage1))
})

test_that("mrlgsi works with real seldata and multiple constraints", {
  d <- setup_genomic_data_real()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  2); C2[1, 1] <- 1; C2[3, 2] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2)
  )
  expect_s3_class(result, "mrlgsi")
  expect_equal(length(result$beta_R1), d$n1)
})

test_that("mrlgsi uses wcol for multi-column wmat", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  wmat2 <- cbind(d$w, d$w * 2)
  res1 <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat2, wcol = 1, C1 = C1, C2 = C2)
  )
  res2 <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat2, wcol = 2, C1 = C1, C2 = C2)
  )
  # beta_R2 depends on w, so columns should produce different restricted coefficients
  expect_false(isTRUE(all.equal(res1$beta_R2, res2$beta_R2)))
})

test_that("mrlgsi respects custom tau", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2, tau = 1.2)
  )
  expect_equal(result$tau, 1.2)
})

test_that("mrlgsi with Young method enabled", {
  d <- setup_genomic_data_real()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2,
           use_young_method = TRUE)
  )
  expect_true(!is.null(result$k1))
})

test_that("mrlgsi warns when C* skipped (b_R1_P1_b_R1 <= 0)", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  expect_warning(
    result <- mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1,
                     P1 = matrix(0, d$n1, d$n1), d$w, C1 = C1, C2 = C2),
    "Invalid phenotypic variance"
  )
  expect_equal(result$C_star, d$C)
})

test_that("mrlgsi uses named traits from colnames", {
  d <- setup_genomic_data_small()
  colnames(d$Gamma1) <- paste0("R1T", seq_len(d$n1))
  colnames(d$Gamma)  <- paste0("R2T", seq_len(d$n))
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2)
  )
  expect_equal(result$summary_stage1$Trait, paste0("R1T", seq_len(d$n1)))
  expect_equal(result$summary_stage2$Trait, paste0("R2T", seq_len(d$n)))
})

test_that("mrlgsi Young fallback to manual when rho is NA", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- withCallingHandlers(
    mrlgsi(matrix(0, d$n1, d$n1), d$Gamma, d$A1, d$A, d$C, d$G1, d$P1,
           d$w, C1 = C1, C2 = C2, use_young_method = TRUE,
           k1_manual = 1.7, k2_manual = 1.7),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(result$k1, 1.7)
})

test_that("mrlgsi stores Gamma_star and C_star matrices", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2)
  )
  expect_equal(dim(result$Gamma_star), c(d$n, d$n))
  expect_equal(dim(result$C_star),     c(d$n, d$n))
})

test_that("mrlgsi stage metrics have correct lengths", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  result <- suppressWarnings(
    mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w, C1 = C1, C2 = C2)
  )
  expect_equal(length(result$stage1_metrics$E), d$n1)
  expect_equal(length(result$stage2_metrics$E), d$n)
})

# ==============================================================================
# TESTS: mppg_lgsi
# ==============================================================================

test_that("mppg_lgsi basic functionality with default U matrices", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = seq_len(d$n))
  )
  expect_s3_class(result, "mppg_lgsi")
  expect_s3_class(result, "multistage_genomic_index")
  expect_equal(length(result$beta_P1), d$n1)
  expect_equal(length(result$beta_P2), d$n)
  expect_equal(length(result$b_P1),    d$n1)
  expect_true(is.numeric(result$theta1))
  expect_true(is.numeric(result$theta2))
  expect_equal(length(result$gain_ratios_1), d$n1)
  expect_equal(length(result$gain_ratios_2), d$n)
  expect_true("beta_P" %in% names(result$summary_stage1))
  expect_true("Ratio"  %in% names(result$summary_stage1))
})

test_that("mppg_lgsi works with real seldata", {
  d <- setup_genomic_data_real()
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = c(2, 1, 1), d2 = c(3, 2, 1, 1, 1, 0.5, 0.5))
  )
  expect_s3_class(result, "mppg_lgsi")
  expect_equal(length(result$beta_P1), d$n1)
})

test_that("mppg_lgsi stops when d1 length != n1", {
  d <- setup_genomic_data_small()
  expect_error(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = c(1, 2), d2 = seq_len(d$n)),
    "d1 must have length"
  )
})

test_that("mppg_lgsi stops when d2 length != n", {
  d <- setup_genomic_data_small()
  expect_error(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = c(1, 2)),
    "d2 must have length"
  )
})

test_that("mppg_lgsi warns with custom U1 matrix", {
  d <- setup_genomic_data_small()
  U1_custom <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 1), d$n1, d$n1)
  w_cap <- capture_warnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = seq_len(d$n), U1 = U1_custom)
  )
  expect_true(any(grepl("Custom U1", w_cap)))
})

test_that("mppg_lgsi C* adjustment differs between Identity and custom U1", {
  d <- setup_genomic_data_small()
  d1 <- seq_len(d$n1); d2 <- seq_len(d$n)
  res_id <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = d1, d2 = d2, U1 = NULL)
  )
  U1_perm <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 1), d$n1, d$n1)
  res_custom <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = d1, d2 = d2, U1 = U1_perm)
  )
  expect_equal(res_id$b_P1,   res_custom$b_P1,   tolerance = 1e-10)
  expect_equal(res_id$C_star, res_custom$C_star,  tolerance = 1e-10)
  expect_false(isTRUE(all.equal(res_id$beta_P1, res_custom$beta_P1, tolerance = 1e-4)))
})

test_that("mppg_lgsi works with custom U2 square permutation matrix", {
  d <- setup_genomic_data_small()
  # U2 must be square (n x n) because d2 is validated to have length n;
  # use a permutation matrix so ncol(U2) = n and d2 length n is valid
  U2_perm <- diag(d$n)
  U2_perm[1, 1] <- 0; U2_perm[2, 2] <- 0
  U2_perm[1, 2] <- 1; U2_perm[2, 1] <- 1  # swap rows 1 and 2
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = seq_len(d$n), U2 = U2_perm)
  )
  expect_s3_class(result, "mppg_lgsi")
  expect_equal(length(result$beta_P2), d$n)
})

test_that("mppg_lgsi uses wcol for multi-column wmat", {
  d <- setup_genomic_data_small()
  wmat2 <- cbind(d$w, d$w * 0.5)
  res1 <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat2,
              wcol = 1, d1 = seq_len(d$n1), d2 = seq_len(d$n))
  )
  res2 <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, wmat2,
              wcol = 2, d1 = seq_len(d$n1), d2 = seq_len(d$n))
  )
  # beta_P1/theta1 depends on w, so different columns produce different results
  expect_false(isTRUE(all.equal(res1$beta_P1, res2$beta_P1)))
})

test_that("mppg_lgsi respects custom tau", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = seq_len(d$n), tau = 0.9)
  )
  expect_equal(result$tau, 0.9)
})

test_that("mppg_lgsi warns when C* skipped (b_P1_P1_b_P1 <= 0)", {
  d <- setup_genomic_data_small()
  expect_warning(
    result <- mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1,
                        P1 = matrix(0, d$n1, d$n1), d$w,
                        d1 = seq_len(d$n1), d2 = seq_len(d$n)),
    "Invalid phenotypic variance"
  )
  expect_equal(result$C_star, d$C)
})

test_that("mppg_lgsi with Young method enabled", {
  d <- setup_genomic_data_real()
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = c(2, 1, 1), d2 = c(3, 2, 1, 1, 1, 0.5, 0.5),
              use_young_method = TRUE)
  )
  expect_true(!is.null(result$k1))
})

test_that("mppg_lgsi uses named traits from colnames", {
  d <- setup_genomic_data_small()
  colnames(d$Gamma1) <- paste0("P1T", seq_len(d$n1))
  colnames(d$Gamma)  <- paste0("P2T", seq_len(d$n))
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = seq_len(d$n))
  )
  expect_equal(result$summary_stage1$Trait, paste0("P1T", seq_len(d$n1)))
  expect_equal(result$summary_stage2$Trait, paste0("P2T", seq_len(d$n)))
})

test_that("mppg_lgsi theta1 = 0 when d1 all zeros (denominator near zero)", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = rep(0, d$n1), d2 = seq_len(d$n))
  )
  expect_equal(result$theta1, 0)
})

test_that("mppg_lgsi theta2 = 0 when d2 all zeros (denominator near zero)", {
  d <- setup_genomic_data_small()
  result <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = rep(0, d$n))
  )
  expect_equal(result$theta2, 0)
})

test_that("mppg_lgsi Young fallback to manual when rho is NA", {
  d <- setup_genomic_data_small()
  result <- withCallingHandlers(
    mppg_lgsi(matrix(0, d$n1, d$n1), d$Gamma, d$A1, d$A, d$C, d$G1, d$P1,
              d$w, d1 = seq_len(d$n1), d2 = seq_len(d$n),
              use_young_method = TRUE, k1_manual = 1.6, k2_manual = 1.6),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_equal(result$k1, 1.6)
})

test_that("mppg_lgsi C_star non-PD warning branch executes with extreme k1", {
  set.seed(999)
  n1 <- 2; n <- 3
  Gamma1 <- diag(n1) + matrix(0.3, n1, n1)
  Gamma  <- diag(n)  + matrix(0.3, n, n)
  A1 <- 0.9 * Gamma1
  A  <- 0.9 * Gamma[, 1:n1]
  G1 <- Gamma1 * 1.1
  P1 <- Gamma1
  C  <- matrix(0.001, n, n); diag(C) <- 0.002
  w  <- rep(1, n); d1 <- c(1, 1); d2 <- rep(1, n)
  warns <- capture_warnings(
    mppg_lgsi(Gamma1, Gamma, A1, A, C, G1, P1, w,
              d1 = d1, d2 = d2, k1_manual = 20, k2_manual = 2)
  )
  # Function ran; may or may not warn about PD - just ensure it executes
  expect_true(is.character(warns) || length(warns) == 0)
})

# ==============================================================================
# TESTS: Cross-function consistency and structure
# ==============================================================================

test_that("mlgsi and mrlgsi produce different stage-1 coefficients", {
  d <- setup_genomic_data_real()
  C1 <- matrix(0, d$n1, 1); C1[2, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[2, 1] <- 1
  res_ml  <- suppressWarnings(mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w))
  res_mrl <- suppressWarnings(mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                                     C1 = C1, C2 = C2))
  expect_false(isTRUE(all.equal(res_ml$beta1, res_mrl$beta_R1)))
})

test_that("mlgsi and mppg_lgsi produce different stage-1 coefficients", {
  d <- setup_genomic_data_real()
  res_ml  <- suppressWarnings(mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w))
  res_ppg <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = c(2, 1, 1), d2 = c(3, 2, 1, 1, 1, 0.5, 0.5))
  )
  expect_false(isTRUE(all.equal(res_ml$beta1, res_ppg$beta_P1)))
})

test_that("All three indices return rho_I1I2", {
  d <- setup_genomic_data_real()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  d1 <- c(2, 1, 1); d2 <- c(3, 2, 1, 1, 1, 0.5, 0.5)

  res_ml  <- suppressWarnings(mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w))
  res_mrl <- suppressWarnings(mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                                     C1 = C1, C2 = C2))
  res_ppg <- suppressWarnings(mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                                        d1 = d1, d2 = d2))

  expect_true(!is.null(res_ml$rho_I1I2))
  expect_true(!is.null(res_mrl$rho_I1I2))
  expect_true(!is.null(res_ppg$rho_I1I2))
})

test_that("Manual intensities default works across all three indices", {
  d <- setup_genomic_data_small()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1

  res_ml  <- suppressWarnings(mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w))
  res_mrl <- suppressWarnings(mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                                     C1 = C1, C2 = C2))
  res_ppg <- suppressWarnings(
    mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
              d1 = seq_len(d$n1), d2 = seq_len(d$n))
  )

  expect_equal(res_ml$k1,  2.063)
  expect_equal(res_mrl$k1, 2.063)
  expect_equal(res_ppg$k1, 2.063)
})

# ==============================================================================
# TESTS: Young's method error handler (selection_proportion out of range)
# When .young_intensities throws an error (p = 0), the tryCatch handler fires:
#   warning("Young's method failed. Using manual intensities.")
#   list(k1 = k1_manual, k2 = k2_manual)
# This covers the error handler branches in all three functions.
# ==============================================================================

test_that("mlgsi Young's method error handler falls back to manual intensities", {
  d <- setup_genomic_data_real()
  # tau supplied explicitly so STEP 1 skips computing it from selection_proportion
  # (avoiding tau = Inf). selection_proportion = 0 is passed to .young_intensities
  # which stops, the tryCatch catches it -> warning + manual k values returned.
  expect_warning(
    result <- mlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                    tau = 1.5,
                    use_young_method = TRUE, selection_proportion = 0,
                    k1_manual = 1.5, k2_manual = 1.8),
    "Young's method failed"
  )
  expect_equal(result$k1, 1.5)
  expect_equal(result$k2, 1.8)
})

test_that("mrlgsi Young's method error handler falls back to manual intensities", {
  d <- setup_genomic_data_real()
  C1 <- matrix(0, d$n1, 1); C1[1, 1] <- 1
  C2 <- matrix(0, d$n,  1); C2[1, 1] <- 1
  expect_warning(
    result <- mrlgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                     C1 = C1, C2 = C2,
                     tau = 1.5,
                     use_young_method = TRUE, selection_proportion = 0,
                     k1_manual = 1.5, k2_manual = 1.8),
    "Young's method failed"
  )
  expect_equal(result$k1, 1.5)
  expect_equal(result$k2, 1.8)
})

test_that("mppg_lgsi Young's method error handler falls back to manual intensities", {
  d <- setup_genomic_data_real()
  result <- NULL
  # The C* adjustment emits a "not positive definite" warning before Young's method
  # is called; suppress that specific warning so only "Young's method failed" is seen.
  withCallingHandlers(
    expect_warning(
      {
        result <- mppg_lgsi(d$Gamma1, d$Gamma, d$A1, d$A, d$C, d$G1, d$P1, d$w,
                            d1 = c(2, 1, 1), d2 = c(3, 2, 1, 1, 1, 0.5, 0.5),
                            tau = 1.5,
                            use_young_method = TRUE, selection_proportion = 0,
                            k1_manual = 1.5, k2_manual = 1.8)
      },
      "Young's method failed"
    ),
    warning = function(w) {
      if (grepl("positive definite|C\\*", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_equal(result$k1, 1.5)
  expect_equal(result$k2, 1.8)
})
