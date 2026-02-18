# ==============================================================================
# Tests for eigen_indices.R — Chapter 7 Eigen Selection Index Methods
#
# Functions covered: esim(), resim(), ppg_esim()
# Helpers tested indirectly: .esim_solve_sym_multi(), .eigen_index_metrics(),
#                             .leading_eigenvector()
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared fixtures  (2-trait synthetic + real seldata matrices)
# ------------------------------------------------------------------------------

# Small 2×2 system with known spectral properties:
#   P^{-1}G = [5 1; 1 4] [10 2; 2 8]^{-1}  — exact eigenvalues computable by hand
P_2  <- matrix(c(10, 2, 2, 8),  nrow = 2)
G_2  <- matrix(c(5,  1, 1, 4),  nrow = 2)

# 3×3 system used for restriction tests
P_3  <- matrix(c(10, 2, 1,
                  2, 8, 1,
                  1, 1, 6), nrow = 3, byrow = TRUE)
G_3  <- matrix(c(5, 1, 0.5,
                  1, 4, 0.5,
                  0.5, 0.5, 3), nrow = 3, byrow = TRUE)
colnames(P_3) <- colnames(G_3) <- c("T1", "T2", "T3")

# ==============================================================================
# =====  ESIM  =================================================================
# ==============================================================================

test_that("esim: output structure is correct", {
  r <- esim(P_2, G_2)

  expect_s3_class(r, "esim")
  expect_s3_class(r, "eigen_index")

  expected_names <- c("summary", "b", "Delta_G", "sigma_I",
                       "hI2", "rHI", "lambda2", "implied_w",
                       "all_eigenvalues", "selection_intensity",
                       "trait_names", "pmat", "gmat")
  expect_true(all(expected_names %in% names(r)))

  # summary is a data frame
  expect_s3_class(r$summary, "data.frame")
  expect_equal(nrow(r$summary), 1L)        # default n_indices = 1

  # b and Delta_G are numeric vectors of correct length
  expect_true(is.numeric(r$b))
  expect_false(is.matrix(r$b))
  expect_equal(length(r$b), 2L)
  expect_equal(length(r$Delta_G), 2L)

  # scalar outputs are single finite numbers
  for (field in c("sigma_I", "hI2", "rHI", "lambda2")) {
    expect_true(is.numeric(r[[field]]) && length(r[[field]]) == 1L,
                label = paste("scalar:", field))
    expect_true(is.finite(r[[field]]), label = paste("finite:", field))
  }
})

test_that("esim: leading eigenvalue equals index heritability", {
  r <- esim(P_2, G_2)

  # The leading eigenvalue of P^{-1}G IS h^2_I by definition
  P_inv_G <- solve(P_2) %*% G_2
  ev       <- eigen(P_inv_G)
  lambda2_ref <- max(Re(ev$values))

  expect_equal(r$lambda2, lambda2_ref, tolerance = 1e-6)
  expect_equal(r$hI2,     lambda2_ref, tolerance = 1e-6)
  expect_equal(r$rHI,     sqrt(lambda2_ref), tolerance = 1e-6)
})

test_that("esim: b_E is the eigenvector of P^{-1}G for the leading eigenvalue", {
  r <- esim(P_2, G_2)

  P_inv_G <- solve(P_2) %*% G_2
  # (P^{-1}G) b = lambda^2 b  =>  residual P^{-1}G b - lambda^2 b ≈ 0
  residual <- as.numeric(P_inv_G %*% r$b - r$lambda2 * r$b)
  expect_true(all(abs(residual) < 1e-8),
              label = "b_E satisfies (P^{-1}G - lambda^2 I) b = 0")
})

test_that("esim: metric formulas are internally consistent", {
  r <- esim(P_3, G_3)
  b <- as.numeric(r$b)

  # sigma_I = sqrt(b'Pb)
  bPb <- as.numeric(t(b) %*% P_3 %*% b)
  expect_equal(r$sigma_I, sqrt(bPb), tolerance = 1e-8)

  # hI2 = b'Gb / b'Pb  (also == lambda2)
  bGb <- as.numeric(t(b) %*% G_3 %*% b)
  expect_equal(r$hI2, bGb / bPb, tolerance = 1e-8)

  # rHI = sqrt(hI2)
  expect_equal(r$rHI, sqrt(r$hI2), tolerance = 1e-8)

  # Delta_G_scalar = k_I * sigma_I
  k_I <- r$selection_intensity
  expect_equal(r$summary$Delta_G[1], k_I * r$sigma_I, tolerance = 1e-6)

  # Delta_G vector = (k_I / sigma_I) * G b
  DG_expected <- as.numeric(k_I * (G_3 %*% b) / r$sigma_I)
  expect_equal(as.numeric(r$Delta_G), DG_expected, tolerance = 1e-8)
})

test_that("esim: n_indices > 1 returns multiple rows in summary", {
  r <- esim(P_3, G_3, n_indices = 3L)

  # P_3 is 3×3, so at most 3 positive eigenvalues
  expect_true(nrow(r$summary) >= 1L)
  expect_true(nrow(r$summary) <= 3L)

  # Eigenvalues in summary should be in descending order
  lambdas <- r$summary$lambda2
  expect_true(all(diff(lambdas) <= 0 + 1e-10))
})

test_that("esim: returns correct trait names from column names", {
  r <- esim(P_3, G_3)
  expect_equal(r$trait_names, c("T1", "T2", "T3"))
  expect_equal(names(r$b),      c("T1", "T2", "T3"))
  expect_equal(names(r$Delta_G), c("T1", "T2", "T3"))
})

test_that("esim: fallback trait names when matrices have no colnames", {
  P_nc <- P_2; colnames(P_nc) <- NULL
  G_nc <- G_2; colnames(G_nc) <- NULL
  r    <- esim(P_nc, G_nc)
  expect_equal(r$trait_names, c("Trait_1", "Trait_2"))
})

test_that("esim: input validation catches bad inputs", {
  expect_error(esim(P_2, G_3),              regexp = "same dimensions")
  expect_error(esim(matrix(1), matrix(1)),  regexp = "At least 2 traits")
  non_sym <- P_2; non_sym[1, 2] <- 999
  expect_error(esim(non_sym, G_2),          regexp = "symmetric")
})

test_that("esim: works on real seldata matrices", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  r <- esim(pmat, gmat)
  expect_s3_class(r, "esim")
  expect_equal(length(r$b), 7L)
  expect_true(r$hI2 > 0 && r$hI2 < 1)
  expect_true(r$rHI > 0 && r$rHI < 1)
  expect_true(r$sigma_I > 0)
})

test_that("esim: print and summary return invisibly without error", {
  r <- esim(P_3, G_3)
  capture.output(expect_invisible(print(r)))
  capture.output(expect_invisible(summary(r)))
})

test_that("esim: results are deterministic across repeated calls", {
  r1 <- esim(P_3, G_3)
  r2 <- esim(P_3, G_3)
  expect_equal(r1$b,       r2$b)
  expect_equal(r1$lambda2, r2$lambda2)
  expect_equal(r1$Delta_G, r2$Delta_G)
})

test_that("esim: custom selection_intensity scales sigma_I correctly", {
  k1 <- 2.063
  k2 <- 1.755    # 20 % selection
  r1 <- esim(P_3, G_3, selection_intensity = k1)
  r2 <- esim(P_3, G_3, selection_intensity = k2)

  # b and lambda2 are independent of k_I (purely from eigenproblem)
  expect_equal(r1$b,       r2$b,       tolerance = 1e-10)
  expect_equal(r1$lambda2, r2$lambda2, tolerance = 1e-10)

  # Delta_G (scalar) scales with k_I
  expect_equal(r1$summary$Delta_G[1] / r2$summary$Delta_G[1], k1 / k2,
               tolerance = 1e-6)
})

# ==============================================================================
# =====  RESIM  ================================================================
# ==============================================================================

test_that("resim: output structure is correct", {
  r <- resim(P_3, G_3, restricted_traits = 1L)

  expect_s3_class(r, "resim")
  expect_s3_class(r, "eigen_index")

  expected_names <- c("summary", "b", "Delta_G", "sigma_I",
                       "hI2", "rHI", "lambda2", "K", "Q_R",
                       "U_mat", "restricted_traits", "implied_w",
                       "selection_intensity", "trait_names", "pmat", "gmat")
  expect_true(all(expected_names %in% names(r)))

  expect_s3_class(r$summary, "data.frame")
  expect_equal(nrow(r$summary), 1L)
  expect_equal(length(r$b), 3L)
  expect_equal(length(r$Delta_G), 3L)

  for (field in c("sigma_I", "hI2", "rHI", "lambda2")) {
    expect_true(is.finite(r[[field]]), label = paste("finite:", field))
  }
})

test_that("resim: single restricted trait has near-zero genetic gain", {
  r <- resim(P_3, G_3, restricted_traits = 1L)
  expect_true(abs(r$Delta_G["T1"]) < 1e-6,
              label = "Restricted trait T1 has Delta_G ≈ 0")
})

test_that("resim: two restricted traits both have near-zero gain", {
  r <- resim(P_3, G_3, restricted_traits = c(1L, 2L))
  expect_true(abs(r$Delta_G["T1"]) < 1e-6)
  expect_true(abs(r$Delta_G["T2"]) < 1e-6)
})

test_that("resim: unrestricted traits are free to respond", {
  r <- resim(P_3, G_3, restricted_traits = 1L)
  # T2 and T3 are unrestricted — their gains should generally be non-zero
  expect_true(abs(r$Delta_G["T2"]) > 1e-8 || abs(r$Delta_G["T3"]) > 1e-8)
})

test_that("resim: restriction satisfied via formula U' G b ≈ 0", {
  r <- resim(P_3, G_3, restricted_traits = c(1L, 3L))
  U  <- r$U_mat
  b  <- as.numeric(r$b)

  # The formal RESIM constraint is U' C b = 0  (C = gmat)
  constraint <- as.numeric(t(U) %*% G_3 %*% b)
  expect_true(all(abs(constraint) < 1e-6),
              label = "U' G b = 0 for all restricted traits")
})

test_that("resim: accepts custom U_mat and matches restricted_traits version", {
  r1 <- resim(P_3, G_3, restricted_traits = c(1L, 2L))
  U  <- diag(3)[, c(1L, 2L), drop = FALSE]
  r2 <- resim(P_3, G_3, U_mat = U)

  expect_equal(r1$b,       r2$b,       tolerance = 1e-8)
  expect_equal(r1$lambda2, r2$lambda2, tolerance = 1e-8)
})

test_that("resim: K is a projection matrix (K^2 = K)", {
  r <- resim(P_3, G_3, restricted_traits = 1L)
  K  <- r$K
  expect_equal(K %*% K, K, tolerance = 1e-8)
})

test_that("resim: lambda2 < esim lambda2 (restriction reduces heritability)", {
  re  <- esim(P_3, G_3)
  rr  <- resim(P_3, G_3, restricted_traits = 1L)
  expect_true(rr$lambda2 <= re$lambda2 + 1e-10,
              label = "Restriction cannot improve maximum index heritability")
})

test_that("resim: input validation catches bad inputs", {
  expect_error(resim(P_3, G_3),
               regexp = "restricted_traits.*U_mat")
  expect_error(resim(P_3, G_3, restricted_traits = c(1L, 2L, 3L)),
               regexp = "less than number of traits")
  expect_error(resim(P_3, G_3, restricted_traits = 5L),
               regexp = "valid trait indices")
})

test_that("resim: works on real seldata matrices", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  r <- resim(pmat, gmat, restricted_traits = c(1L, 3L))
  expect_s3_class(r, "resim")
  expect_true(abs(r$Delta_G[1]) < 1e-6)   # trait 1 restricted
  expect_true(abs(r$Delta_G[3]) < 1e-6)   # trait 3 restricted
})

test_that("resim: print and summary return invisibly without error", {
  r <- resim(P_3, G_3, restricted_traits = 1L)
  capture.output(expect_invisible(print(r)))
  capture.output(expect_invisible(summary(r)))
})

# ==============================================================================
# =====  PPG-ESIM  =============================================================
# ==============================================================================

test_that("ppg_esim: output structure is correct", {
  d <- c(1, 1)
  r <- ppg_esim(P_2, G_2, d)

  expect_s3_class(r, "ppg_esim")
  expect_s3_class(r, "eigen_index")

  expected_names <- c("summary", "beta", "b", "Delta_G", "sigma_I",
                       "hI2", "rHI", "lambda2", "F_mat", "K_P",
                       "Q_P", "D_M", "desired_gains",
                       "selection_intensity", "trait_names", "pmat", "gmat")
  expect_true(all(expected_names %in% names(r)))

  expect_s3_class(r$summary, "data.frame")
  expect_equal(length(r$beta), 2L)
  expect_equal(length(r$b),    2L)
  expect_equal(length(r$Delta_G), 2L)

  for (field in c("sigma_I", "hI2", "rHI", "lambda2")) {
    expect_true(is.finite(r[[field]]), label = paste("finite:", field))
  }
})

test_that("ppg_esim: Mallard matrix D_M is orthogonal to d", {
  d <- c(2, 1, -1)
  r <- ppg_esim(P_3, G_3, d)

  d_unit <- d / sqrt(sum(d^2))
  # Every column of D_M must be orthogonal to d_unit
  inner_products <- as.numeric(t(r$D_M) %*% d_unit)
  expect_true(all(abs(inner_products) < 1e-10),
              label = "D_M columns are orthogonal to d")

  # D_M columns should be orthonormal among themselves
  DtD <- t(r$D_M) %*% r$D_M
  expect_equal(DtD, diag(ncol(r$D_M)), tolerance = 1e-10)
})

test_that("ppg_esim: K_P is a rank-1 projection matrix", {
  d <- c(1, 2, 1)
  r <- ppg_esim(P_3, G_3, d)

  K <- r$K_P
  # Idempotent: K^2 = K
  expect_equal(K %*% K, K, tolerance = 1e-8)

  # Rank 1
  rank_K <- sum(eigen(K, only.values = TRUE)$values > 0.5)
  expect_equal(rank_K, 1L)
})

test_that("ppg_esim: gains are exactly proportional to d (equal d)", {
  d <- rep(1, 3)
  r <- ppg_esim(P_3, G_3, d)

  dg  <- as.numeric(r$Delta_G)
  ratios <- dg / d
  cv   <- sd(ratios) / abs(mean(ratios))
  expect_true(cv < 1e-8,
              label = paste("CV of gain ratios:", cv))
})

test_that("ppg_esim: gains are exactly proportional to d (asymmetric d)", {
  d <- c(3, 1, 2)
  r <- ppg_esim(P_3, G_3, d)

  dg     <- as.numeric(r$Delta_G)
  ratios <- dg / d
  cv     <- sd(ratios) / abs(mean(ratios))
  expect_true(cv < 1e-8,
              label = paste("CV of gain ratios:", cv))
})

test_that("ppg_esim: gains are exactly proportional to d (mixed sign d)", {
  d <- c(1, -1, 1)
  r <- ppg_esim(P_3, G_3, d)

  dg     <- as.numeric(r$Delta_G)
  ratios <- dg / d
  cv     <- sd(ratios) / abs(mean(ratios))
  expect_true(cv < 1e-8,
              label = paste("CV of gain ratios:", cv))
})

test_that("ppg_esim: gains are exactly proportional on real seldata (7 traits)", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  for (dtest in list(rep(1, 7), c(2, 1, 1, 1, 1, 1, 1), c(1, 1, -1, 1, -1, 1, 1))) {
    r      <- ppg_esim(pmat, gmat, dtest)
    dg     <- as.numeric(r$Delta_G)
    ratios <- dg / dtest
    cv     <- sd(ratios) / abs(mean(ratios))
    expect_true(cv < 1e-6,
                label = paste("CV:", cv, "for d:", paste(dtest, collapse = ",")))
  }
})

test_that("ppg_esim: phi (mean ratio) is positive — gains align with d", {
  d1 <- c(1, 1, 1)
  d2 <- c(2, 1, -1)

  for (d in list(d1, d2)) {
    r   <- ppg_esim(P_3, G_3, d)
    dg  <- as.numeric(r$Delta_G)
    phi <- mean(dg / d)
    expect_true(phi > 0,
                label = paste("phi > 0 for d:", paste(d, collapse = ",")))
  }
})

test_that("ppg_esim: beta = F b  and F is scalar ±I", {
  d <- c(1, 2, 1)
  r <- ppg_esim(P_3, G_3, d)

  # beta = F_mat %*% b
  beta_reconstructed <- as.numeric(r$F_mat %*% r$b)
  expect_equal(as.numeric(r$beta), beta_reconstructed, tolerance = 1e-10)

  # F_mat must be a scalar multiple of identity (±I)
  diag_vals <- diag(r$F_mat)
  expect_true(length(unique(diag_vals)) == 1L,
              label = "F_mat is scalar ±I")
  expect_true(abs(abs(diag_vals[1]) - 1) < 1e-10,
              label = "diagonal of F_mat is ±1")
  # Off-diagonal elements are zero
  off_diag <- r$F_mat - diag(diag_vals)
  expect_true(all(abs(off_diag) < 1e-10))
})

test_that("ppg_esim: metric formulas are internally consistent", {
  d <- c(2, 1, 1)
  r <- ppg_esim(P_3, G_3, d)

  beta <- as.numeric(r$beta)

  bPb     <- as.numeric(t(beta) %*% P_3 %*% beta)
  bGb     <- as.numeric(t(beta) %*% G_3 %*% beta)
  k_I     <- r$selection_intensity

  expect_equal(r$sigma_I, sqrt(bPb), tolerance = 1e-8)
  expect_equal(r$hI2,     bGb / bPb, tolerance = 1e-8)
  expect_equal(r$rHI,     sqrt(r$hI2), tolerance = 1e-8)

  DG_expected <- as.numeric(k_I * (G_3 %*% beta) / r$sigma_I)
  expect_equal(as.numeric(r$Delta_G), DG_expected, tolerance = 1e-8)
})

test_that("ppg_esim: D_M has correct dimensions (t x t-1)", {
  d <- c(1, 2, 3)
  r <- ppg_esim(P_3, G_3, d)
  expect_equal(dim(r$D_M), c(3L, 2L))
})

test_that("ppg_esim: scaling d does not change the gain direction", {
  d1 <- c(1, 2, 1)
  d2 <- d1 * 5          # scaled version
  r1 <- ppg_esim(P_3, G_3, d1)
  r2 <- ppg_esim(P_3, G_3, d2)

  # Gains should be proportional (same direction, possibly different scale)
  ratio1 <- as.numeric(r1$Delta_G) / as.numeric(r2$Delta_G)
  cv     <- sd(ratio1) / abs(mean(ratio1))
  expect_true(cv < 1e-8,
              label = "Scaling d does not change the gain ratio pattern")
})

test_that("ppg_esim: lambda2 <= esim lambda2 (PPG constraint reduces heritability)", {
  d  <- c(2, 1, 1)
  re <- esim(P_3, G_3)
  rp <- ppg_esim(P_3, G_3, d)
  expect_true(rp$lambda2 <= re$lambda2 + 1e-10,
              label = "PPG constraint cannot improve free-optimum heritability")
})

test_that("ppg_esim: input validation catches bad inputs", {
  expect_error(ppg_esim(P_3, G_3, c(1, 1)),
               regexp = "same length")
  expect_error(ppg_esim(P_3, G_3, c(0, 0, 0)),
               regexp = "non-zero")
  non_sym <- P_3; non_sym[1, 2] <- 999
  expect_error(ppg_esim(non_sym, G_3, c(1, 1, 1)),
               regexp = "symmetric")
})

test_that("ppg_esim: print and summary return invisibly without error", {
  d <- c(2, 1, 1)
  r <- ppg_esim(P_3, G_3, d)
  capture.output(expect_invisible(print(r)))
  capture.output(expect_invisible(summary(r)))
})

test_that("ppg_esim: results are deterministic across repeated calls", {
  d <- c(1, 2, 1)
  r1 <- ppg_esim(P_3, G_3, d)
  r2 <- ppg_esim(P_3, G_3, d)
  expect_equal(r1$beta,    r2$beta)
  expect_equal(r1$lambda2, r2$lambda2)
  expect_equal(r1$Delta_G, r2$Delta_G)
})

# ==============================================================================
# =====  Cross-method consistency  =============================================
# ==============================================================================

test_that("esim >= resim >= ppg_esim in heritability (hierarchy of constraints)", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  re  <- esim(pmat, gmat)
  rr  <- resim(pmat, gmat, restricted_traits = 1L)
  rp  <- ppg_esim(pmat, gmat, rep(1, 7))

  # Free ESIM achieves the maximum; restrictions can only reduce it
  expect_true(re$lambda2 >= rr$lambda2 - 1e-8)
  expect_true(re$lambda2 >= rp$lambda2 - 1e-8)
})

test_that("esim without restrictions equals resim with zero restrictions (edge case)", {
  # RESIM with restriction on 0 traits is ill-formed — verify proper error
  # (not a valid call — at least 1 restriction required)
  expect_error(resim(P_3, G_3, restricted_traits = integer(0)),
               regexp = "restricted_traits")
})
