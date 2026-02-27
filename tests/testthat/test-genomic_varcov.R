# ==============================================================================
# Tests for Genomic Variance-Covariance Functions (R/genomic.R)
# ==============================================================================

# ==============================================================================
# TEST: genomic_varcov() - Genomic Variance-Covariance Matrix (Γ)
# ==============================================================================

test_that("genomic_varcov returns correct dimensions", {
  set.seed(123)
  n_genotypes <- 100
  n_traits <- 5
  gebv_mat <- matrix(rnorm(n_genotypes * n_traits),
    nrow = n_genotypes, ncol = n_traits
  )
  colnames(gebv_mat) <- paste0("Trait", 1:n_traits)

  Gamma <- genomic_varcov(gebv_mat)

  expect_equal(nrow(Gamma), n_traits)
  expect_equal(ncol(Gamma), n_traits)
  expect_equal(dim(Gamma), c(n_traits, n_traits))
})

test_that("genomic_varcov returns symmetric matrix", {
  set.seed(456)
  gebv_mat <- matrix(rnorm(500), nrow = 100, ncol = 5)

  Gamma <- genomic_varcov(gebv_mat)

  expect_true(isSymmetric(Gamma))
  expect_true(isSymmetric(unname(Gamma)))
})

test_that("genomic_varcov preserves trait names", {
  set.seed(789)
  gebv_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)
  trait_names <- c("Yield", "Height", "Maturity", "Quality", "Biomass")
  colnames(gebv_mat) <- trait_names

  Gamma <- genomic_varcov(gebv_mat)

  expect_equal(colnames(Gamma), trait_names)
  expect_equal(rownames(Gamma), trait_names)
})

test_that("genomic_varcov works with different correlation methods", {
  set.seed(111)
  gebv_mat <- matrix(rnorm(300), nrow = 60, ncol = 5)

  Gamma_pearson <- genomic_varcov(gebv_mat, method = "pearson")
  Gamma_kendall <- genomic_varcov(gebv_mat, method = "kendall")
  Gamma_spearman <- genomic_varcov(gebv_mat, method = "spearman")

  expect_equal(dim(Gamma_pearson), c(5, 5))
  expect_equal(dim(Gamma_kendall), c(5, 5))
  expect_equal(dim(Gamma_spearman), c(5, 5))

  # Methods should give different results
  expect_false(identical(Gamma_pearson, Gamma_kendall))
})

test_that("genomic_varcov handles missing values with complete.obs", {
  set.seed(222)
  gebv_mat <- matrix(rnorm(300), nrow = 60, ncol = 5)
  gebv_mat[c(1, 5, 10), 1] <- NA
  gebv_mat[c(3, 7), 3] <- NA

  Gamma <- genomic_varcov(gebv_mat, use = "complete.obs")

  expect_equal(dim(Gamma), c(5, 5))
  expect_true(all(is.finite(Gamma)))
  expect_true(isSymmetric(Gamma))
})

test_that("genomic_varcov warns with pairwise.complete.obs and missing values", {
  set.seed(333)
  gebv_mat <- matrix(rnorm(300), nrow = 60, ncol = 5)
  gebv_mat[c(1, 5, 10), 1] <- NA

  expect_warning(
    Gamma <- genomic_varcov(gebv_mat, use = "pairwise.complete.obs"),
    "Missing values detected"
  )

  expect_equal(dim(Gamma), c(5, 5))
})

test_that("genomic_varcov errors with everything and missing values", {
  set.seed(444)
  gebv_mat <- matrix(rnorm(300), nrow = 60, ncol = 5)
  gebv_mat[1, 1] <- NA

  expect_error(
    genomic_varcov(gebv_mat, use = "everything"),
    "Missing values detected"
  )
})

test_that("genomic_varcov errors with non-numeric data", {
  gebv_mat <- matrix(letters[1:20], nrow = 4, ncol = 5)

  expect_error(
    genomic_varcov(gebv_mat),
    "numeric"
  )
})

test_that("genomic_varcov errors with too few observations", {
  gebv_mat <- matrix(rnorm(5), nrow = 1, ncol = 5)

  expect_error(
    genomic_varcov(gebv_mat),
    "at least 2 observations"
  )
})

test_that("genomic_varcov returns all finite values", {
  set.seed(555)
  gebv_mat <- matrix(rnorm(500), nrow = 100, ncol = 5)

  Gamma <- genomic_varcov(gebv_mat)

  expect_true(all(is.finite(Gamma)))
  expect_false(any(is.na(Gamma)))
  expect_false(any(is.infinite(Gamma)))
})

# ==============================================================================
# TEST: phenomic_genomic_varcov() - Phenomic-Genomic Covariance Matrix (Φ)
# ==============================================================================

test_that("phenomic_genomic_varcov returns correct dimensions from raw data", {
  set.seed(666)
  n_genotypes <- 100
  n_traits <- 7

  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )
  gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
    nrow = n_genotypes, ncol = n_traits
  )

  Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)

  # Should be 2*n_traits x 2*n_traits
  expect_equal(dim(Phi), c(2 * n_traits, 2 * n_traits))
  expect_equal(nrow(Phi), 14)
  expect_equal(ncol(Phi), 14)
})

test_that("phenomic_genomic_varcov returns symmetric matrix", {
  set.seed(777)
  phen_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)
  gebv_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)

  Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)

  expect_true(isSymmetric(unname(Phi), tol = 1e-10))
})

test_that("phenomic_genomic_varcov works with pre-computed matrices", {
  set.seed(888)
  n_traits <- 5

  P <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)
  P <- (P + t(P)) / 2 # Make symmetric

  Gamma <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)
  Gamma <- (Gamma + t(Gamma)) / 2 # Make symmetric

  P_yg <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)

  Phi <- phenomic_genomic_varcov(P = P, Gamma = Gamma, P_yg = P_yg)

  expect_equal(dim(Phi), c(2 * n_traits, 2 * n_traits))
  expect_true(isSymmetric(unname(Phi), tol = 1e-10))
})

test_that("phenomic_genomic_varcov adds dimension names", {
  set.seed(999)
  phen_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)
  gebv_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)
  colnames(phen_mat) <- paste0("T", 1:5)
  colnames(gebv_mat) <- paste0("T", 1:5)

  Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)

  expect_false(is.null(colnames(Phi)))
  expect_false(is.null(rownames(Phi)))
  expect_true(any(grepl("_phen", colnames(Phi))))
  expect_true(any(grepl("_gebv", colnames(Phi))))
})

test_that("phenomic_genomic_varcov errors without necessary inputs", {
  expect_error(
    phenomic_genomic_varcov(),
    "Must provide either"
  )
})

test_that("phenomic_genomic_varcov errors with dimension mismatch", {
  phen_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)
  gebv_mat <- matrix(rnorm(300), nrow = 60, ncol = 5) # Different rows

  expect_error(
    phenomic_genomic_varcov(phen_mat, gebv_mat),
    "same number of rows"
  )

  gebv_mat2 <- matrix(rnorm(400), nrow = 80, ncol = 5)
  phen_mat2 <- matrix(rnorm(320), nrow = 80, ncol = 4) # Different cols

  expect_error(
    phenomic_genomic_varcov(phen_mat2, gebv_mat2),
    "same number of columns"
  )
})

test_that("phenomic_genomic_varcov errors with non-symmetric P or Gamma", {
  n_traits <- 5

  P <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)
  # P is not symmetric

  Gamma <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)
  Gamma <- (Gamma + t(Gamma)) / 2

  P_yg <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)

  expect_error(
    phenomic_genomic_varcov(P = P, Gamma = Gamma, P_yg = P_yg),
    "symmetric"
  )
})

test_that("phenomic_genomic_varcov works with different covariance methods", {
  set.seed(1010)
  phen_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)
  gebv_mat <- matrix(rnorm(400), nrow = 80, ncol = 5)

  Phi_pearson <- phenomic_genomic_varcov(phen_mat, gebv_mat, method = "pearson")
  Phi_kendall <- phenomic_genomic_varcov(phen_mat, gebv_mat, method = "kendall")

  expect_equal(dim(Phi_pearson), dim(Phi_kendall))
  expect_false(identical(Phi_pearson, Phi_kendall))
})

# ==============================================================================
# TEST: genetic_genomic_varcov() - Genetic-Genomic Covariance Matrix (A)
# ==============================================================================

test_that("genetic_genomic_varcov returns correct dimensions", {
  set.seed(1111)
  n_traits <- 7

  gmat <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)
  gmat <- (gmat + t(gmat)) / 2 # Make symmetric

  A <- genetic_genomic_varcov(gmat)

  # Should be 2*n_traits x 2*n_traits
  expect_equal(dim(A), c(2 * n_traits, 2 * n_traits))
  expect_equal(nrow(A), 14)
})

test_that("genetic_genomic_varcov returns symmetric matrix", {
  set.seed(1212)
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  A <- genetic_genomic_varcov(gmat)

  expect_true(isSymmetric(unname(A), tol = 1e-10))
})

test_that("genetic_genomic_varcov works with reliability parameter", {
  set.seed(1313)
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Single reliability value
  A1 <- genetic_genomic_varcov(gmat, reliability = 0.8)
  expect_equal(nrow(A1), 2 * ncol(gmat))

  # Vector of reliabilities
  rel_vec <- seq(0.6, 0.9, length.out = ncol(gmat))
  A2 <- genetic_genomic_varcov(gmat, reliability = rel_vec)
  expect_equal(nrow(A2), 2 * ncol(gmat))
})

test_that("genetic_genomic_varcov works with Gamma parameter", {
  set.seed(1414)
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  Gamma <- gmat * 0.8 # Simulate genomic variance

  A <- genetic_genomic_varcov(gmat, Gamma = Gamma)

  expect_equal(nrow(A), 2 * ncol(gmat))
  expect_true(isSymmetric(unname(A), tol = 1e-10))
})

test_that("genetic_genomic_varcov works with C_gebv_g parameter", {
  set.seed(1515)
  n_traits <- 5
  gmat <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)
  gmat <- (gmat + t(gmat)) / 2

  C_gebv_g <- matrix(rnorm(n_traits^2), nrow = n_traits, ncol = n_traits)

  A <- genetic_genomic_varcov(gmat, C_gebv_g = C_gebv_g)

  expect_equal(dim(A), c(2 * n_traits, 2 * n_traits))
})

test_that("genetic_genomic_varcov square parameter works", {
  set.seed(1616)
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])


  A_square <- genetic_genomic_varcov(gmat, square = TRUE)
  expect_equal(nrow(A_square), 2 * ncol(gmat))
  expect_equal(ncol(A_square), 2 * ncol(gmat))


  A_rect <- genetic_genomic_varcov(gmat, square = FALSE)
  expect_equal(nrow(A_rect), 2 * ncol(gmat))
  expect_equal(ncol(A_rect), ncol(gmat))
})

test_that("genetic_genomic_varcov errors with invalid reliability", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])


  expect_error(
    genetic_genomic_varcov(gmat, reliability = 1.5),
    "reliability|0 and 1"
  )


  expect_error(
    genetic_genomic_varcov(gmat, reliability = -0.2),
    "reliability|0 and 1"
  )

  # Wrong length vector
  expect_error(
    genetic_genomic_varcov(gmat, reliability = c(0.7, 0.8)),
    "reliability|length"
  )
})

test_that("genetic_genomic_varcov errors with dimension mismatch", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  wrong_Gamma <- matrix(rnorm(25), nrow = 5, ncol = 5)
  wrong_Gamma <- (wrong_Gamma + t(wrong_Gamma)) / 2 # Make symmetric

  expect_error(
    genetic_genomic_varcov(gmat, Gamma = wrong_Gamma),
    "Gamma must be"
  )
})

test_that("genetic_genomic_varcov adds dimension names", {
  set.seed(1717)
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  A <- genetic_genomic_varcov(gmat)

  expect_false(is.null(colnames(A)))
  expect_false(is.null(rownames(A)))
})

test_that("genetic_genomic_varcov returns all finite values", {
  set.seed(1818)
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  A <- genetic_genomic_varcov(gmat)

  expect_true(all(is.finite(A)))
  expect_false(any(is.na(A)))
})

# ==============================================================================
# TEST: Integration tests
# ==============================================================================

test_that("genomic covariance functions work together", {
  set.seed(2020)

  # Generate test data
  n_genotypes <- 100
  n_traits <- 7 # Match seldata trait count

  phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
    nrow = n_genotypes, ncol = n_traits
  )
  gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
    nrow = n_genotypes, ncol = n_traits
  )

  # Compute all matrices
  Gamma <- genomic_varcov(gebv_mat)
  Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)

  # Use realvariance-covariance matrices
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  A <- genetic_genomic_varcov(gmat, Gamma = Gamma[seq_len(nrow(gmat)), seq_len(ncol(gmat))]) # 1:nrow(gmat), 1:ncol(gmat)])

  expect_equal(dim(Gamma), c(n_traits, n_traits))
  expect_equal(dim(Phi), c(2 * n_traits, 2 * n_traits))
  expect_equal(nrow(A), 2 * nrow(gmat))
})

# ==============================================================================
# NEW COVERAGE TESTS — targeting previously uncovered lines in R/genomic.R
# ==============================================================================

test_that("phenomic_genomic_varcov parameter validation and warnings (lines 241-275)", {
  n_traits <- 3
  P <- matrix(rnorm(n_traits^2), n_traits, n_traits)
  P <- (P + t(P)) / 2 # symmetric

  Gamma_nonsym <- matrix(rnorm(n_traits^2), n_traits, n_traits)
  # line 241: Gamma not symmetric
  expect_error(
    phenomic_genomic_varcov(P = P, Gamma = Gamma_nonsym, P_yg = P),
    "Gamma must be symmetric"
  )

  Gamma_sym <- (Gamma_nonsym + t(Gamma_nonsym)) / 2
  Gamma_wrong_dim <- matrix(1, 2, 2)
  # line 247: Gamma wrong dimension
  expect_error(
    phenomic_genomic_varcov(P = P, Gamma = Gamma_wrong_dim, P_yg = P),
    "Gamma must be 3 x 3"
  )

  P_yg_wrong_dim <- matrix(1, 2, 2)
  # line 250: P_yg wrong dimension
  expect_error(
    phenomic_genomic_varcov(P = P, Gamma = Gamma_sym, P_yg = P_yg_wrong_dim),
    "P_yg must be 3 x 3"
  )

  # lines 272-275: Phi not symmetric
  # Mock is_symmetric to return FALSE when called on the combined 6x6 Phi matrix
  mock_is_sym <- function(x, ...) {
    if (nrow(x) > n_traits) {
      return(FALSE)
    }
    TRUE
  }

  expect_warning(
    with_mocked_bindings(
      phenomic_genomic_varcov(P = P, Gamma = Gamma_sym, P_yg = P),
      is_symmetric = mock_is_sym,
      .package = "selection.index"
    ),
    "Phi is not symmetric \\(max asymmetry"
  )
})

test_that("genetic_genomic_varcov parameter validation and warnings (lines 366-439)", {
  n_traits <- 3
  gmat_nonsym <- matrix(rnorm(n_traits^2), n_traits, n_traits)

  # line 366: gmat not symmetric
  expect_error(
    genetic_genomic_varcov(gmat = gmat_nonsym),
    "gmat must be symmetric"
  )

  gmat_sym <- (gmat_nonsym + t(gmat_nonsym)) / 2
  Gamma_nonsym <- matrix(rnorm(n_traits^2), n_traits, n_traits)

  # line 376: Gamma not symmetric
  expect_error(
    genetic_genomic_varcov(gmat = gmat_sym, Gamma = Gamma_nonsym),
    "Gamma must be symmetric"
  )

  Gamma_sym <- (Gamma_nonsym + t(Gamma_nonsym)) / 2
  C_gebv_g_wrong <- matrix(1, 2, 2)

  # line 388: C_gebv_g wrong dimension
  expect_error(
    genetic_genomic_varcov(gmat = gmat_sym, Gamma = Gamma_sym, C_gebv_g = C_gebv_g_wrong),
    "C_gebv_g must be 3 x 3"
  )

  # lines 436-439: A is not symmetric
  mock_is_sym <- function(x, ...) {
    if (nrow(x) > n_traits) {
      return(FALSE)
    }
    TRUE
  }

  expect_warning(
    with_mocked_bindings(
      genetic_genomic_varcov(gmat = gmat_sym, Gamma = Gamma_sym),
      is_symmetric = mock_is_sym,
      .package = "selection.index"
    ),
    "A is not symmetric \\(max asymmetry"
  )
})
