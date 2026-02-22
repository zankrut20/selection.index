# test-phenotypic_indices.R
# Comprehensive tests for R/phenotypic_indices.R targeting 95%+ coverage
# Functions: smith_hazel, base_index, lpsi (combinatorial)
#            print.smith_hazel, summary.smith_hazel
#            print.base_index, summary.base_index
#            .index_metrics (internal, covered indirectly)

# ==============================================================================
# TEST DATA SETUP
# ==============================================================================

# Small synthetic data (3 traits) with known structure
setup_phen_data_small <- function(n = 3, seed = 42) {
  set.seed(seed)
  # Construct PD matrices via Cholesky
  L <- matrix(c(
    2, 0, 0,
    1, 2, 0,
    0.5, 1, 1.5
  ), nrow = n, byrow = TRUE)
  P <- t(L) %*% L
  colnames(P) <- rownames(P) <- paste0("t", seq_len(n))

  L2 <- matrix(c(
    1, 0, 0,
    0.6, 1, 0,
    0.2, 0.4, 0.8
  ), nrow = n, byrow = TRUE)
  G <- t(L2) %*% L2
  colnames(G) <- rownames(G) <- paste0("t", seq_len(n))

  w <- as.numeric(seq_len(n))
  list(P = P, G = G, w = w, n = n)
}

# Real data using seldata package dataset
setup_phen_data_real <- function() {
  data("seldata", package = "selection.index", envir = environment())
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  w <- as.numeric(c(10, 8, 6, 4, 2, 1, 1)) # 7-trait weights
  list(
    pmat = pmat, gmat = gmat, w = w, n = 7,
    trait_names = colnames(pmat)
  )
}

# ==============================================================================
# .index_metrics  (covered indirectly through smith_hazel / base_index)
# ==============================================================================

test_that(".index_metrics NA branches fire when G = zero matrix", {
  d <- setup_phen_data_small()
  G_zero <- matrix(0, d$n, d$n)
  colnames(G_zero) <- rownames(G_zero) <- paste0("t", seq_len(d$n))

  # b = P^{-1} * G_zero * w = 0 → bPb = 0 → all metrics NA
  res <- smith_hazel(d$P, G_zero, d$w)
  expect_true(is.na(res$sigma_I))
  expect_true(is.na(res$GA))
  expect_true(is.na(res$PRE))
  expect_true(is.na(res$hI2))
  expect_true(is.na(res$rHI))
  expect_true(all(is.na(res$Delta_G)))
})

test_that(".index_metrics PRE_constant uses GAY when provided", {
  d <- setup_phen_data_small()
  res_no_GAY <- smith_hazel(d$P, d$G, d$w)
  res_with_GAY <- smith_hazel(d$P, d$G, d$w, GAY = 10)
  # Without GAY: PRE = GA * 100; with GAY: PRE = GA * 100/10 = GA * 10
  expect_equal(res_with_GAY$PRE, res_no_GAY$PRE / 10, tolerance = 1e-6)
})

test_that(".index_metrics GA is 0 (not NA) in base_index when G = zero (b = w, bPb > 0)", {
  d <- setup_phen_data_small()
  G_zero <- matrix(0, d$n, d$n)
  colnames(G_zero) <- rownames(G_zero) <- paste0("t", seq_len(d$n))
  # base_index: b = w (not zero), bPb = w'Pw > 0, sigma_I > 0
  # bGw = w' 0 w = 0  => GA = i * 0 / sigma_I = 0  (finite, not NA)
  res <- base_index(d$P, G_zero, d$w, compare_to_lpsi = FALSE)
  expect_equal(res$GA, 0)
  expect_equal(res$PRE, 0)
})

# ==============================================================================
# SMITH-HAZEL INDEX
# ==============================================================================

test_that("smith_hazel returns correct structure with small synthetic data", {
  d <- setup_phen_data_small()
  res <- smith_hazel(d$P, d$G, d$w)

  expect_s3_class(res, "smith_hazel")
  expect_s3_class(res, "selection_index")
  expect_named(res, c(
    "b", "w", "Delta_G", "sigma_I", "GA", "PRE",
    "hI2", "rHI", "selection_intensity", "summary"
  ))
  expect_length(res$b, d$n)
  expect_length(res$w, d$n)
  expect_length(res$Delta_G, d$n)
  expect_true(is.numeric(res$sigma_I))
  expect_true(res$sigma_I > 0)
  expect_true(res$rHI > 0 && res$rHI <= 1)
  expect_true(res$hI2 >= 0)
  expect_true(is.data.frame(res$summary))
})

test_that("smith_hazel works with real seldata", {
  d <- setup_phen_data_real()
  res <- smith_hazel(d$pmat, d$gmat, d$w)

  expect_s3_class(res, "smith_hazel")
  expect_length(res$b, d$n)
  expect_gt(res$GA, 0)
  expect_gt(res$rHI, 0)
  expect_equal(names(res$w), d$trait_names)
})

test_that("smith_hazel handles matrix wmat with wcol selection", {
  d <- setup_phen_data_real()
  wmat2 <- cbind(d$w, rev(d$w))
  res1 <- smith_hazel(d$pmat, d$gmat, wmat2, wcol = 1)
  res2 <- smith_hazel(d$pmat, d$gmat, wmat2, wcol = 2)

  expect_equal(res1$w, setNames(d$w, d$trait_names))
  expect_equal(res2$w, setNames(rev(d$w), d$trait_names))
  expect_false(identical(res1$b, res2$b))
})

test_that("smith_hazel respects custom selection_intensity", {
  d <- setup_phen_data_small()
  res_default <- smith_hazel(d$P, d$G, d$w)
  res_custom <- smith_hazel(d$P, d$G, d$w, selection_intensity = 1.755)

  expect_equal(res_custom$selection_intensity, 1.755)
  expect_equal(res_custom$b, res_default$b) # b doesn't change with intensity
  expect_equal(res_custom$sigma_I, res_default$sigma_I)
  expect_equal(res_custom$GA / res_default$GA, 1.755 / 2.063, tolerance = 1e-6)
})

test_that("smith_hazel handles GAY for PRE calculation", {
  d <- setup_phen_data_small()
  res <- smith_hazel(d$P, d$G, d$w, GAY = 5)
  expect_true(is.finite(res$PRE))
  expect_equal(res$PRE, res$GA * 100 / 5, tolerance = 1e-6)
})

test_that("smith_hazel stops when pmat is not square", {
  expect_error(
    smith_hazel(matrix(1:6, 2, 3), diag(2), c(1, 2)),
    "square"
  )
})

test_that("smith_hazel stops when gmat is not square", {
  expect_error(
    smith_hazel(diag(2), matrix(1:6, 2, 3), c(1, 2)),
    "square"
  )
})

test_that("smith_hazel stops when pmat and gmat have different dimensions", {
  expect_error(
    smith_hazel(diag(3), diag(2), c(1, 2)),
    "same dimensions"
  )
})

test_that("smith_hazel stops when wmat has wrong number of rows", {
  d <- setup_phen_data_small()
  expect_error(
    smith_hazel(d$P, d$G, wmat = c(1, 2)), # 2 instead of 3
    "Number of rows in wmat"
  )
})

test_that("smith_hazel stops when wcol is out of range", {
  d <- setup_phen_data_small()
  expect_error(
    smith_hazel(d$P, d$G, cbind(d$w, d$w), wcol = 5),
    "wcol must be between"
  )
  expect_error(
    smith_hazel(d$P, d$G, d$w, wcol = 0),
    "wcol must be between"
  )
})

test_that("smith_hazel stops on non-finite economic weights", {
  d <- setup_phen_data_small()
  w_bad <- c(1, Inf, 3)
  expect_error(
    smith_hazel(d$P, d$G, w_bad),
    "finite"
  )
  w_nan <- c(1, NaN, 3)
  expect_error(
    smith_hazel(d$P, d$G, w_nan),
    "finite"
  )
})

# ==============================================================================
# BASE INDEX
# ==============================================================================

test_that("base_index returns correct structure with small synthetic data", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w)

  expect_s3_class(res, "base_index")
  expect_s3_class(res, "selection_index")
  expect_named(res, c(
    "b", "w", "Delta_G", "sigma_I", "GA", "PRE",
    "hI2", "rHI", "selection_intensity", "summary",
    "lpsi_comparison"
  ))
  # In Base Index b = w
  expect_equal(res$b, round(d$w, 4))
  expect_true(res$sigma_I > 0)
  expect_false(is.null(res$lpsi_comparison))
})

test_that("base_index works with real seldata", {
  d <- setup_phen_data_real()
  res <- base_index(d$pmat, d$gmat, d$w)

  expect_s3_class(res, "base_index")
  expect_length(res$b, d$n)
  expect_gt(res$GA, 0)
  expect_false(is.null(res$lpsi_comparison))
  expect_true(is.numeric(res$lpsi_comparison$efficiency_ratio))
})

test_that("base_index with compare_to_lpsi = FALSE omits comparison", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w, compare_to_lpsi = FALSE)

  expect_null(res$lpsi_comparison)
})

test_that("base_index handles matrix wmat with wcol", {
  d <- setup_phen_data_small()
  wmat2 <- cbind(d$w, rev(d$w))
  res1 <- base_index(d$P, d$G, wmat2, wcol = 1)
  res2 <- base_index(d$P, d$G, wmat2, wcol = 2)

  expect_equal(res1$b, round(d$w, 4))
  expect_equal(res2$b, round(rev(d$w), 4))
})

test_that("base_index handles GAY for PRE calculation", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w, GAY = 5)
  expect_equal(res$PRE, res$GA * 100 / 5, tolerance = 1e-6)
})

test_that("base_index LPSI comparison tryCatch error handler produces warning on failure", {
  d <- setup_phen_data_small()
  # cpp_symmetric_solve is called by .solve_sym_multi inside the tryCatch block.
  # Mock it to throw so the error handler fires, emitting the warning and leaving
  # lpsi_comparison NULL.
  testthat::with_mocked_bindings(
    cpp_symmetric_solve = function(...) stop("mocked solver failure"),
    .package = "selection.index",
    code = {
      expect_warning(
        res <- base_index(d$P, d$G, d$w, compare_to_lpsi = TRUE),
        "Could not calculate LPSI comparison"
      )
      expect_null(res$lpsi_comparison)
    }
  )
})

test_that("base_index efficiency_ratio is NA when GA_lpsi <= 0", {
  d <- setup_phen_data_small()
  G_zero <- matrix(0, d$n, d$n)
  colnames(G_zero) <- rownames(G_zero) <- paste0("t", seq_len(d$n))
  # LPSI GA = 0 → efficiency ratio NA
  res <- base_index(d$P, G_zero, d$w, compare_to_lpsi = TRUE)
  expect_true(is.na(res$lpsi_comparison$efficiency_ratio))
})

test_that("base_index stops on non-square matrices", {
  expect_error(
    base_index(matrix(1:6, 2, 3), diag(2), c(1, 2)),
    "square"
  )
})

test_that("base_index stops when pmat and gmat differ in dimension", {
  expect_error(
    base_index(diag(3), diag(2), c(1, 2)),
    "same dimensions"
  )
})

test_that("base_index stops on wrong wmat row count", {
  d <- setup_phen_data_small()
  expect_error(
    base_index(d$P, d$G, c(1, 2)),
    "Number of rows in wmat"
  )
})

test_that("base_index stops on out-of-range wcol", {
  d <- setup_phen_data_small()
  expect_error(
    base_index(d$P, d$G, d$w, wcol = 0),
    "wcol must be between"
  )
})

test_that("base_index stops on non-finite weights", {
  d <- setup_phen_data_small()
  expect_error(
    base_index(d$P, d$G, c(1, NA, 3)),
    "finite"
  )
})

# ==============================================================================
# COMBINATORIAL LPSI
# ==============================================================================

test_that("lpsi returns a data frame with expected columns", {
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  res <- lpsi(ncomb = 3, pmat = d$pmat, gmat = d$gmat, wmat = wmat, wcol = 1)

  expect_true(is.data.frame(res))
  expect_true("ID" %in% names(res))
  expect_true("GA" %in% names(res))
  expect_true("PRE" %in% names(res))
  expect_true("rHI" %in% names(res))
  expect_true("hI2" %in% names(res))
  expect_true("Rank" %in% names(res))
  expect_equal(nrow(res), choose(d$n, 3))
})

test_that("lpsi with GAY produces correct PRE_constant", {
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  res_no_GAY <- lpsi(3, d$pmat, d$gmat, wmat, wcol = 1)
  res_with_GAY <- lpsi(3, d$pmat, d$gmat, wmat, wcol = 1, GAY = 5)

  # With GAY: PRE = GA * 100/5 = GA * 20; without GAY: PRE = GA * 100
  expect_equal(res_with_GAY$PRE, res_no_GAY$PRE / 5, tolerance = 1e-4)
})

test_that("lpsi excluding_trait as numeric vector filters correctly", {
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  # Exclude trait 1 and 2 from combinations
  res <- lpsi(3, d$pmat, d$gmat, wmat, wcol = 1, excluding_trait = c(1, 2))

  # No combo should contain index 1 or 2
  combos <- strsplit(res$ID, ", ")
  has_excluded <- sapply(combos, function(ids) any(ids %in% c("1", "2")))
  expect_false(any(has_excluded))
  expect_lt(nrow(res), choose(d$n, 3))
})

test_that("lpsi excluding_trait as character vector filters correctly", {
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  # Exclude by trait name
  res <- lpsi(3, d$pmat, d$gmat, wmat,
    wcol = 1,
    excluding_trait = c("sypp", "dtf")
  )

  expect_true(is.data.frame(res))
  expect_lt(nrow(res), choose(d$n, 3))
})

test_that("lpsi excluding_trait character warns when no trait names match", {
  d <- setup_phen_data_small()
  wmat <- matrix(d$w, ncol = 1)
  expect_warning(
    res <- lpsi(2, d$P, d$G, wmat, excluding_trait = c("nonexistent")),
    "None of the specified trait names found"
  )
  # No filtering → all combinations returned
  expect_equal(nrow(res), choose(d$n, 2))
})

test_that("lpsi excluding_trait character stops when pmat has no colnames", {
  d <- setup_phen_data_small()
  P_nonames <- d$P
  colnames(P_nonames) <- NULL
  wmat <- matrix(d$w, ncol = 1)
  expect_error(
    lpsi(2, P_nonames, d$G, wmat, excluding_trait = c("t1")),
    "pmat must have column names"
  )
})

test_that("lpsi excluding_trait as data.frame filters correctly", {
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  # Create a data frame whose column names match trait names to exclude
  excl_df <- as.data.frame(matrix(0, nrow = 1, ncol = 2))
  colnames(excl_df) <- c("sypp", "dtf")

  res <- lpsi(3, d$pmat, d$gmat, wmat, wcol = 1, excluding_trait = excl_df)
  expect_true(is.data.frame(res))
  expect_lt(nrow(res), choose(d$n, 3))
})

test_that("lpsi excluding_trait data.frame warns when no column names match", {
  d <- setup_phen_data_small()
  wmat <- matrix(d$w, ncol = 1)
  excl_df <- data.frame(zzz = 1)
  expect_warning(
    res <- lpsi(2, d$P, d$G, wmat, excluding_trait = excl_df),
    "None of the column names from excluding_trait found"
  )
  expect_equal(nrow(res), choose(d$n, 2))
})

test_that("lpsi excluding_trait data.frame stops when pmat has no colnames", {
  d <- setup_phen_data_small()
  P_nonames <- d$P
  colnames(P_nonames) <- NULL
  wmat <- matrix(d$w, ncol = 1)
  excl_df <- data.frame(t1 = 1)
  expect_error(
    lpsi(2, P_nonames, d$G, wmat, excluding_trait = excl_df),
    "pmat must have column names"
  )
})

test_that("lpsi stops when excluding_trait is invalid type", {
  d <- setup_phen_data_small()
  wmat <- matrix(d$w, ncol = 1)
  expect_error(
    lpsi(2, d$P, d$G, wmat, excluding_trait = TRUE),
    "must be a numeric vector"
  )
})

test_that("lpsi returns empty data frame when all combinations excluded", {
  d <- setup_phen_data_small()
  wmat <- matrix(d$w, ncol = 1)
  # Exclude all traits → no valid combination of size 2 exists
  res <- lpsi(2, d$P, d$G, wmat, excluding_trait = c(1, 2, 3))

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 0)
  expect_true("ID" %in% names(res))
})

test_that("lpsi Rank column uses min ties", {
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  res <- lpsi(2, d$pmat, d$gmat, wmat, wcol = 1)

  expect_true(all(res$Rank >= 1))
  expect_equal(min(res$Rank), 1)
})

test_that("lpsi rHI and hI2 are 0 when G = zero matrix (bPb = 0 else branches)", {
  # With G = zero: Gw_full = 0, b = P^{-1}*0 = 0, bPb = 0, wGw_full = 0
  # → hI2 else branch fires (= 0) and rHI else branch fires (= 0, line 646)
  d <- setup_phen_data_small()
  G_zero <- matrix(0, d$n, d$n)
  colnames(G_zero) <- rownames(G_zero) <- paste0("t", seq_len(d$n))
  wmat <- matrix(d$w, ncol = 1)
  res <- suppressWarnings(lpsi(2, d$P, G_zero, wmat))

  expect_true(all(res$rHI == 0))
  expect_true(all(res$hI2 == 0))
})

test_that("lpsi excluding_trait matrix without colnames triggers stop (via data.frame-like path)", {
  # Note: is.numeric(matrix(1,...)) is TRUE in R, so numeric matrices fall into
  # Case 1 (numeric index). The 'excluding_trait data must have column names' stop
  # is reachable only via a non-numeric, non-character matrix; in practice this is
  # effectively dead code given how R dispatches is.numeric/is.character on matrices.
  # Verify numeric matrix IS treated as numeric index (no error):
  d <- setup_phen_data_small()
  wmat <- matrix(d$w, ncol = 1)
  excl_mat <- matrix(1L, nrow = 1, ncol = 1) # numeric matrix: treated as index c(1)
  res <- lpsi(2, d$P, d$G, wmat, excluding_trait = excl_mat)
  # Trait 1 excluded from 2-trait combos: only (2,3) remains
  expect_equal(nrow(res), 1L)
})

# ==============================================================================
# PRINT.SMITH_HAZEL
# ==============================================================================

test_that("print.smith_hazel produces expected output with named traits", {
  d <- setup_phen_data_real()
  res <- smith_hazel(d$pmat, d$gmat, d$w)

  out <- capture.output(print(res))
  expect_true(any(grepl("SMITH-HAZEL", out)))
  expect_true(any(grepl("Genetic Advance", out)))
  expect_true(any(grepl("Accuracy", out)))
  expect_true(any(grepl("EXPECTED GENETIC RESPONSE", out)))
})

test_that("print.smith_hazel uses Trait_N when trait names are missing", {
  # pmat without colnames → names(x$w) is NULL → fallback to Trait_1, Trait_2, ...
  P <- diag(3)
  G <- diag(3)
  w <- c(1, 2, 3)
  res <- smith_hazel(P, G, w)
  # names(w) will be NULL since P has no colnames

  out <- capture.output(print(res))
  expect_true(any(grepl("Trait_1", out)))
})

test_that("print.smith_hazel skips PRE line when PRE is NA", {
  d <- setup_phen_data_small()
  G_zero <- matrix(0, d$n, d$n)
  colnames(G_zero) <- rownames(G_zero) <- paste0("t", seq_len(d$n))
  res <- smith_hazel(d$P, G_zero, d$w)
  # GA and PRE are NA → the PRE line must NOT appear
  out <- capture.output(print(res))
  expect_false(any(grepl("Relative Efficiency", out)))
})

test_that("print.smith_hazel returns invisible(x)", {
  d <- setup_phen_data_small()
  res <- smith_hazel(d$P, d$G, d$w)
  capture.output(ret <- withVisible(print(res)))
  expect_false(ret$visible)
  expect_identical(ret$value, res)
})

# ==============================================================================
# SUMMARY.SMITH_HAZEL
# ==============================================================================

test_that("summary.smith_hazel prints additional statistics", {
  d <- setup_phen_data_real()
  res <- smith_hazel(d$pmat, d$gmat, d$w)

  out <- capture.output(summary(res))
  expect_true(any(grepl("ADDITIONAL STATISTICS", out)))
  expect_true(any(grepl("Economic Weights", out)))
  expect_true(any(grepl("Expected Genetic Gains", out)))
  expect_true(any(grepl("Index Coefficients", out)))
})

test_that("summary.smith_hazel returns invisible(object)", {
  d <- setup_phen_data_small()
  res <- smith_hazel(d$P, d$G, d$w)
  capture.output(ret <- withVisible(summary(res)))
  expect_false(ret$visible)
  expect_identical(ret$value, res)
})

# ==============================================================================
# PRINT.BASE_INDEX
# ==============================================================================

test_that("print.base_index produces expected output with named traits", {
  d <- setup_phen_data_real()
  res <- base_index(d$pmat, d$gmat, d$w)

  out <- capture.output(print(res))
  expect_true(any(grepl("BASE INDEX", out)))
  expect_true(any(grepl("Genetic Advance", out)))
  expect_true(any(grepl("COMPARISON WITH OPTIMAL LPSI", out)))
})

test_that("print.base_index uses Trait_N when trait names are missing", {
  P <- diag(3)
  G <- diag(3)
  w <- c(1, 2, 3)
  res <- base_index(P, G, w, compare_to_lpsi = FALSE)

  out <- capture.output(print(res))
  expect_true(any(grepl("Trait_1", out)))
})

test_that("print.base_index skips LPSI section when compare_to_lpsi = FALSE", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w, compare_to_lpsi = FALSE)

  out <- capture.output(print(res))
  expect_false(any(grepl("COMPARISON WITH OPTIMAL LPSI", out)))
  # PRE itself is still printed (GA is finite when b=w with PD P)
  # but the LPSI comparison block is absent
  expect_false(any(grepl("Efficiency Ratio", out)))
})

test_that("print.base_index skips PRE line when PRE is NA via smith_hazel with G_zero", {
  # In base_index, b=w and P is PD → bPb > 0 → PRE is never NA.
  # The PRE=NA branch in print is reachable only via smith_hazel with G_zero.
  d <- setup_phen_data_small()
  G_zero <- matrix(0, d$n, d$n)
  colnames(G_zero) <- rownames(G_zero) <- paste0("t", seq_len(d$n))
  res <- smith_hazel(d$P, G_zero, d$w)
  expect_true(is.na(res$PRE))

  out <- capture.output(print(res))
  expect_false(any(grepl("Relative Efficiency", out)))
})

test_that("print.base_index shows efficiency_ratio >= 0.95 message", {
  # Identity matrices → b_lpsi = b_base = w → efficiency_ratio = 1.0
  P <- diag(3)
  G <- diag(3)
  w <- c(1, 2, 3)
  res <- base_index(P, G, w, compare_to_lpsi = TRUE)

  expect_gte(res$lpsi_comparison$efficiency_ratio, 0.95)
  out <- capture.output(print(res))
  expect_true(any(grepl(">=95%", out)))
})

test_that("print.base_index shows efficiency_ratio < 0.9 message", {
  # Contrived: use a mock object with efficiency_ratio = 0.85
  mock_res <- structure(
    list(
      b = c(1, 2, 3),
      w = c(1, 2, 3),
      Delta_G = c(0.5, 1.0, 1.5),
      sigma_I = 2.0,
      GA = 0.85,
      PRE = 85,
      hI2 = 0.7,
      rHI = 0.84,
      selection_intensity = 2.063,
      summary = data.frame(),
      lpsi_comparison = list(
        b_lpsi = c(0.5, 1.5, 2.5),
        GA_lpsi = 1.0,
        PRE_lpsi = 100,
        hI2_lpsi = 0.8,
        rHI_lpsi = 0.89,
        Delta_G_lpsi = c(0.6, 1.1, 1.6),
        efficiency_ratio = 0.85
      )
    ),
    class = c("base_index", "selection_index", "list")
  )
  out <- capture.output(print(mock_res))
  expect_true(any(grepl("<90%", out)))
})

test_that("print.base_index returns invisible(x)", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w)
  capture.output(ret <- withVisible(print(res)))
  expect_false(ret$visible)
  expect_identical(ret$value, res)
})

# ==============================================================================
# SUMMARY.BASE_INDEX
# ==============================================================================

test_that("summary.base_index prints additional details with comparison", {
  d <- setup_phen_data_real()
  res <- base_index(d$pmat, d$gmat, d$w)

  out <- capture.output(summary(res))
  expect_true(any(grepl("ADDITIONAL DETAILS", out)))
  expect_true(any(grepl("LPSI vs Base Index", out)))
  expect_true(any(grepl("Response correlation", out)))
})

test_that("summary.base_index prints without lpsi_comparison when disabled", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w, compare_to_lpsi = FALSE)

  out <- capture.output(summary(res))
  expect_true(any(grepl("ADDITIONAL DETAILS", out)))
  expect_false(any(grepl("LPSI vs Base Index", out)))
})

test_that("summary.base_index shows low-correlation warning when cor < 0.8", {
  # Mock object: Delta_G and Delta_G_lpsi are negatively correlated → cor = -1
  mock_res <- structure(
    list(
      b = c(1, 2, 3),
      w = setNames(c(1, 2, 3), c("t1", "t2", "t3")),
      Delta_G = setNames(c(1, -1, 1), c("t1", "t2", "t3")),
      sigma_I = 1,
      GA = 0.5,
      PRE = 50,
      hI2 = 0.6,
      rHI = 0.77,
      selection_intensity = 2.063,
      summary = data.frame(),
      lpsi_comparison = list(
        b_lpsi = c(2, 2, 2),
        GA_lpsi = 0.6,
        PRE_lpsi = 60,
        hI2_lpsi = 0.7,
        rHI_lpsi = 0.84,
        Delta_G_lpsi = setNames(c(-1, 1, -1), c("t1", "t2", "t3")),
        efficiency_ratio = 0.833
      )
    ),
    class = c("base_index", "selection_index", "list")
  )
  out <- capture.output(summary(mock_res))
  expect_true(any(grepl("Low correlation", out)))
})

test_that("summary.base_index returns invisible(object)", {
  d <- setup_phen_data_small()
  res <- base_index(d$P, d$G, d$w)
  capture.output(ret <- withVisible(summary(res)))
  expect_false(ret$visible)
  expect_identical(ret$value, res)
})

# ==============================================================================
# CROSS-FUNCTION CONSISTENCY
# ==============================================================================

test_that("smith_hazel GA >= base_index GA (LPSI is optimal)", {
  d <- setup_phen_data_real()
  res_sh <- smith_hazel(d$pmat, d$gmat, d$w)
  res_bi <- base_index(d$pmat, d$gmat, d$w, compare_to_lpsi = FALSE)

  # LPSI maximises GA → smith_hazel GA should be >= base_index GA
  expect_gte(res_sh$GA, res_bi$GA - 1e-8)
})

test_that("smith_hazel b matches base_index lpsi_comparison b_lpsi", {
  d <- setup_phen_data_small()
  res_sh <- smith_hazel(d$P, d$G, d$w)
  res_bi <- base_index(d$P, d$G, d$w, compare_to_lpsi = TRUE)

  # Both solve for the same LPSI coefficients (res_sh$b is rounded to 4 dp)
  expect_equal(res_sh$b, res_bi$lpsi_comparison$b_lpsi, tolerance = 1e-3)
})

test_that("lpsi(ncomb=n) top PRE matches smith_hazel PRE", {
  # When using all n traits, top-ranked lpsi combo PRE should match smith_hazel PRE
  d <- setup_phen_data_real()
  data("weight", package = "selection.index", envir = environment())
  wmat <- weight_mat(weight)
  res_lpsi <- lpsi(d$n, d$pmat, d$gmat, wmat, wcol = 1)
  res_sh <- smith_hazel(d$pmat, d$gmat, wmat[, 1], GAY = NULL)

  top_PRE <- res_lpsi$PRE[res_lpsi$Rank == 1]
  # PRE from lpsi uses PRE_constant = 100; smith_hazel also uses 100 when GAY is NULL
  expect_equal(top_PRE, round(res_sh$PRE, 4), tolerance = 0.01)
})

# ==============================================================================
# NEW COVERAGE TESTS — targeting previously uncovered lines
# ==============================================================================

test_that("smith_hazel stops when b coefficients are not finite (line 222)", {
  d <- setup_phen_data_small()

  # Mock cpp_symmetric_solve to return NAs to simulate poorly conditioned matrices
  testthat::with_mocked_bindings(
    cpp_symmetric_solve = function(A, B) {
      rep(NA_real_, length(d$w))
    },
    .package = "selection.index",
    code = {
      expect_error(
        smith_hazel(d$P, d$G, d$w),
        "Failed to compute index coefficients. Check matrix conditioning."
      )
    }
  )
})
