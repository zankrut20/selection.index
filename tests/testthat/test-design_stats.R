# test-design_stats.R
# Comprehensive tests for R/design_stats.R targeting 95%+ coverage

# ==============================================================================
# SHARED TEST DATA HELPERS
# ==============================================================================

# Standard 3x3 RCBD (3 genotypes x 3 reps = 9 obs)
rcbd_data <- function() {
  list(
    trait1  = c(10, 12, 11, 13, 14, 12, 15, 16, 14),
    trait2  = c(15, 16, 14, 18, 17, 15, 19, 20, 18),
    gen_idx = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
    rep_idx = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)
  )
}

# Standard 3x3 LSD (3 treatments, 3 rows, 3 columns)
lsd_data <- function() {
  list(
    trait1  = c(10, 12, 11, 13, 14, 12, 15, 16, 14),
    trait2  = c(15, 16, 14, 18, 17, 15, 19, 20, 18),
    gen_idx = c(1L, 2L, 3L, 2L, 3L, 1L, 3L, 1L, 2L),
    rep_idx = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L),
    col_idx = c(1L, 2L, 3L, 2L, 3L, 1L, 3L, 1L, 2L)
  )
}

# SPD: 3 reps x 2 main plots x 2 sub-plots = 12 obs
spd_data <- function() {
  list(
    trait1   = c(10, 12, 11, 13, 14, 12, 15, 16, 14, 13, 12, 11),
    trait2   = c(15, 16, 14, 18, 17, 15, 19, 20, 18, 17, 16, 15),
    gen_idx  = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
    rep_idx  = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L),
    main_idx = c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L)
  )
}

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

test_that("design_stats stops when LSD called without columns", {
  d <- rcbd_data()
  expect_error(
    selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
      design_type = "LSD"
    ),
    "columns"
  )
})

test_that("design_stats stops when SPD called without main_plots", {
  d <- rcbd_data()
  expect_error(
    selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
      design_type = "SPD"
    ),
    "main_plots"
  )
})

test_that("design_stats coerces non-numeric (integer) traits to numeric", {
  d <- rcbd_data()
  t1_int <- as.integer(d$trait1)
  t2_int <- as.integer(d$trait2)
  res <- selection.index:::design_stats(t1_int, t2_int, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "all"
  )
  expect_true(all(is.finite(c(res$CF, res$GMP, res$EMP))))
})

# ==============================================================================
# RCBD – all four calc_types
# ==============================================================================

test_that("RCBD sums_of_products returns correct fields and values", {
  d <- rcbd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "sums_of_products"
  )
  expect_true(all(c(
    "CF", "TSP", "GSP", "RSP", "ESP", "DFG", "DFR", "DFE",
    "n_genotypes", "n_replications", "design_type"
  ) %in% names(res)))
  expect_true(all(is.finite(c(res$CF, res$TSP, res$GSP, res$RSP, res$ESP))))
  expect_equal(res$design_type, "RCBD")
  expect_equal(res$DFG, 2L)
  expect_equal(res$DFR, 2L)
  expect_equal(res$DFE, 4L)
  # GMP/EMP must NOT be present
  expect_false("GMP" %in% names(res))

  expect_equal(res$TSP, res$GSP + res$RSP + res$ESP, tolerance = 1e-10)
})

test_that("RCBD mean_products returns correct fields", {
  d <- rcbd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "mean_products"
  )
  expect_true(all(c(
    "GMP", "EMP", "DFG", "DFR", "DFE",
    "n_genotypes", "n_replications", "design_type"
  ) %in% names(res)))
  expect_true(all(is.finite(c(res$GMP, res$EMP))))
  expect_false("TSP" %in% names(res))
})

test_that("RCBD anova_stats returns DFs only (no SPs or MPs)", {
  d <- rcbd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "anova_stats"
  )
  expect_true(all(c("DFG", "DFR", "DFE", "n_genotypes", "n_replications", "CF", "design_type") %in% names(res)))
  expect_equal(res$DFG, 2L)
  expect_equal(res$DFR, 2L)
  expect_equal(res$DFE, 4L)
  expect_false("GSP" %in% names(res))
  expect_false("GMP" %in% names(res))
})

test_that("RCBD all returns all fields with correct identities", {
  d <- rcbd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "all"
  )
  expect_true(all(c(
    "CF", "TSP", "GSP", "RSP", "ESP", "GMP", "EMP",
    "DFG", "DFR", "DFE", "n_genotypes", "n_replications", "design_type"
  ) %in% names(res)))
  expect_equal(res$TSP, res$GSP + res$RSP + res$ESP, tolerance = 1e-10)
  expect_equal(res$GMP, res$GSP / res$DFG, tolerance = 1e-10)
  expect_equal(res$EMP, res$ESP / res$DFE, tolerance = 1e-10)
})

test_that("RCBD variance (trait1 == trait2) gives non-negative sums of squares", {
  d <- rcbd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait1, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "all"
  )
  expect_gte(res$GSP, 0)
  expect_gte(res$RSP, 0)
})

# ==============================================================================
# LSD – all four calc_types
# ==============================================================================

test_that("LSD sums_of_products returns correct fields and identity", {
  d <- lsd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    columns = d$col_idx,
    design_type = "LSD", calc_type = "sums_of_products"
  )
  expect_true(all(c(
    "CF", "TSP", "GSP", "RSP", "CSP", "ESP", "DFG", "DFR", "DFC",
    "DFE", "n_genotypes", "n_rows", "n_columns", "design_type"
  ) %in% names(res)))
  expect_equal(res$design_type, "LSD")
  expect_equal(res$DFG, 2L)
  expect_equal(res$DFC, 2L)
  expect_equal(res$DFE, 2L)
  expect_equal(res$TSP, res$GSP + res$RSP + res$CSP + res$ESP, tolerance = 1e-10)
  expect_false("GMP" %in% names(res))
})

test_that("LSD mean_products returns correct fields", {
  d <- lsd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    columns = d$col_idx,
    design_type = "LSD", calc_type = "mean_products"
  )
  expect_true(all(c(
    "GMP", "EMP", "DFG", "DFR", "DFC", "DFE",
    "n_genotypes", "n_rows", "n_columns", "design_type"
  ) %in% names(res)))
  expect_false("TSP" %in% names(res))
  expect_false("CSP" %in% names(res))
})

test_that("LSD anova_stats returns DFs only", {
  d <- lsd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    columns = d$col_idx,
    design_type = "LSD", calc_type = "anova_stats"
  )
  expect_true(all(c(
    "DFG", "DFR", "DFC", "DFE", "n_genotypes", "n_rows",
    "n_columns", "CF", "design_type"
  ) %in% names(res)))
  expect_false("GSP" %in% names(res))
  expect_false("GMP" %in% names(res))
  expect_equal(res$design_type, "LSD")
})

test_that("LSD all returns all components with correct identities", {
  d <- lsd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    columns = d$col_idx,
    design_type = "LSD", calc_type = "all"
  )
  expect_true(all(c(
    "CF", "TSP", "GSP", "RSP", "CSP", "ESP", "GMP", "EMP",
    "DFG", "DFR", "DFC", "DFE", "design_type"
  ) %in% names(res)))
  expect_equal(res$GMP, res$GSP / res$DFG, tolerance = 1e-10)
  expect_equal(res$EMP, res$ESP / res$DFE, tolerance = 1e-10)
})

test_that("LSD variance (trait1 == trait2) gives non-negative sums of squares", {
  d <- lsd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait1, d$gen_idx, d$rep_idx,
    columns = d$col_idx,
    design_type = "LSD", calc_type = "all"
  )
  expect_gte(res$GSP, 0)
  expect_gte(res$RSP, 0)
  expect_gte(res$CSP, 0)
})

# ==============================================================================
# SPD – all four calc_types
# ==============================================================================

test_that("SPD sums_of_products returns correct fields", {
  d <- spd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    main_plots = d$main_idx,
    design_type = "SPD", calc_type = "sums_of_products"
  )
  expected <- c(
    "CF", "TSP", "RSP", "MSP", "GSP", "IMSP", "ESP_MAIN", "ESP",
    "DFR", "DFM", "DFE_MAIN", "DFG", "DFIM", "DFE",
    "n_replications", "n_main_plots", "n_genotypes", "design_type"
  )
  expect_true(all(expected %in% names(res)))
  expect_equal(res$design_type, "SPD")
  expect_false("GMP" %in% names(res))
})

test_that("SPD mean_products returns GMP, EMP, EMP_MAIN", {
  d <- spd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    main_plots = d$main_idx,
    design_type = "SPD", calc_type = "mean_products"
  )
  expect_true(all(c(
    "GMP", "EMP", "EMP_MAIN", "DFR", "DFM", "DFE_MAIN",
    "DFG", "DFIM", "DFE", "n_replications", "n_main_plots",
    "n_genotypes", "design_type"
  ) %in% names(res)))
  expect_false("TSP" %in% names(res))
})

test_that("SPD anova_stats returns DFs and counts only", {
  d <- spd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    main_plots = d$main_idx,
    design_type = "SPD", calc_type = "anova_stats"
  )
  expect_true(all(c(
    "DFR", "DFM", "DFE_MAIN", "DFG", "DFIM", "DFE",
    "n_replications", "n_main_plots", "n_genotypes", "CF", "design_type"
  ) %in% names(res)))
  expect_false("GSP" %in% names(res))
  expect_false("GMP" %in% names(res))
  # r=3, a=2, b=2 → DFR=2, DFM=1, DFE_MAIN=2, DFG=1
  expect_equal(res$DFR, 2L)
  expect_equal(res$DFM, 1L)
  expect_equal(res$DFE_MAIN, 2L)
  expect_equal(res$DFG, 1L)
})

test_that("SPD all returns all components with correct identities", {
  d <- spd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait2, d$gen_idx, d$rep_idx,
    main_plots = d$main_idx,
    design_type = "SPD", calc_type = "all"
  )
  expect_true(all(c(
    "CF", "TSP", "RSP", "MSP", "GSP", "IMSP", "ESP_MAIN", "ESP",
    "GMP", "EMP", "EMP_MAIN", "DFR", "DFM", "DFE_MAIN", "DFG", "DFIM",
    "DFE", "n_replications", "n_main_plots", "n_genotypes", "design_type"
  ) %in% names(res)))
  expect_equal(res$design_type, "SPD")
  expect_equal(res$n_main_plots, 2L)
  expect_equal(res$n_genotypes, 2L)
  expect_equal(res$n_replications, 3L)
  expect_equal(res$GMP, res$GSP / res$DFG, tolerance = 1e-10)
  expect_equal(res$EMP, res$ESP / res$DFE, tolerance = 1e-10)
  expect_equal(res$EMP_MAIN, res$ESP_MAIN / res$DFE_MAIN, tolerance = 1e-10)
})

test_that("SPD variance (trait1 == trait2) gives non-negative main effects", {
  d <- spd_data()
  res <- selection.index:::design_stats(d$trait1, d$trait1, d$gen_idx, d$rep_idx,
    main_plots = d$main_idx,
    design_type = "SPD", calc_type = "all"
  )
  expect_gte(res$GSP, 0)
  expect_gte(res$MSP, 0)
  expect_gte(res$RSP, 0)
})

# ==============================================================================
# REAL DATA INTEGRATION
# ==============================================================================

test_that("RCBD with seldata subset produces correct dimensions and finite values", {
  data("seldata", package = "selection.index", envir = environment())
  idx <- which(seldata[["treat"]] %in% c("G1", "G2", "G3"))
  trait1 <- as.numeric(seldata[["sypp"]][idx])
  trait2 <- as.numeric(seldata[["dtf"]][idx])
  gen_idx <- as.integer(as.factor(seldata[["treat"]][idx]))
  rep_idx <- as.integer(seldata[["rep"]][idx])

  res <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx,
    design_type = "RCBD", calc_type = "all"
  )
  expect_equal(res$n_genotypes, 3L)
  expect_equal(res$n_replications, 3L)
  expect_true(all(is.finite(c(res$CF, res$GMP, res$EMP))))
})

# ==============================================================================
# DESIGN_STATS_API
# ==============================================================================

test_that("design_stats_api RCBD (design_type=1) returns symmetric matrices", {
  d <- rcbd_data()
  data_mat <- cbind(d$trait1, d$trait2)

  res <- selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
    design_type = 1L
  )
  expect_named(res, c(
    "GMS", "EMS", "EMS_MAIN", "DFG", "DFE", "DFE_MAIN",
    "n_rep", "n_gen", "n_main", "MSG", "MSE"
  ), ignore.order = TRUE)
  expect_equal(dim(res$MSG), c(2L, 2L))
  expect_equal(dim(res$MSE), c(2L, 2L))
  expect_equal(res$MSG, t(res$MSG))
  expect_equal(res$MSE, t(res$MSE))
  # Diagonal must match design_stats GMP/EMP for trait1 vs trait1
  ds <- selection.index:::design_stats(d$trait1, d$trait1, d$gen_idx, d$rep_idx,
    design_type = "RCBD", calc_type = "mean_products"
  )
  expect_equal(res$MSG[1, 1], ds$GMP, tolerance = 1e-10)
  expect_equal(res$MSE[1, 1], ds$EMP, tolerance = 1e-10)
  # Non-SPD fields must be NA
  expect_true(is.na(res$DFE_MAIN))
  expect_true(is.na(res$n_main))
})

test_that("design_stats_api LSD (design_type=2) returns correct structure", {
  d <- lsd_data()
  data_mat <- cbind(d$trait1, d$trait2)

  res <- selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
    col_idx = d$col_idx,
    design_type = 2L
  )
  expect_equal(dim(res$MSG), c(2L, 2L))
  expect_equal(dim(res$MSE), c(2L, 2L))
  expect_equal(res$MSG, t(res$MSG))
  expect_equal(res$MSE, t(res$MSE))
  expect_true(all(is.finite(res$GMS)))
  expect_true(all(is.finite(res$EMS)))
  expect_true(is.na(res$DFE_MAIN))
  expect_true(is.na(res$n_main))
})

test_that("design_stats_api SPD (design_type=3) returns main-plot error", {
  d <- spd_data()
  data_mat <- cbind(d$trait1, d$trait2)

  res <- selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
    main_idx = d$main_idx,
    design_type = 3L
  )
  expect_equal(dim(res$MSG), c(2L, 2L))
  expect_equal(dim(res$MSE), c(2L, 2L))
  expect_true(all(is.finite(res$EMS_MAIN)))
  expect_false(is.na(res$DFE_MAIN))
  expect_false(is.na(res$n_main))
  expect_equal(res$n_main, 2L)
  # Diagonal of EMS_MAIN must match EMP_MAIN from design_stats
  ds <- selection.index:::design_stats(d$trait1, d$trait1, d$gen_idx, d$rep_idx,
    main_plots = d$main_idx,
    design_type = "SPD", calc_type = "mean_products"
  )
  expect_equal(res$EMS_MAIN[1], ds$EMP_MAIN, tolerance = 1e-10)
})

test_that("design_stats_api stops on invalid design_type integer", {
  d <- rcbd_data()
  data_mat <- cbind(d$trait1, d$trait2)
  expect_error(
    selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
      design_type = 99L
    ),
    "design_type must be"
  )
})

test_that("design_stats_api single-trait matrix returns 1x1 matrices", {
  d <- rcbd_data()
  data_mat <- matrix(d$trait1, ncol = 1)

  res <- selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
    design_type = 1L
  )
  expect_equal(dim(res$MSG), c(1L, 1L))
  expect_equal(dim(res$MSE), c(1L, 1L))
  expect_length(res$GMS, 1L)
  expect_length(res$EMS, 1L)
})

test_that("design_stats_api 3-trait RCBD produces symmetric 3x3 matrices", {
  d <- rcbd_data()
  data_mat <- cbind(d$trait1, d$trait2, d$trait1 + d$trait2)

  res <- selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
    design_type = 1L
  )
  expect_equal(dim(res$MSG), c(3L, 3L))
  expect_equal(dim(res$MSE), c(3L, 3L))
  expect_equal(res$MSG, t(res$MSG))
  expect_equal(res$MSE, t(res$MSE))
})

test_that("design_stats_api legacy compatibility: MSG diagonal matches GMS", {
  d <- rcbd_data()
  data_mat <- cbind(d$trait1, d$trait2)
  res <- selection.index:::design_stats_api(data_mat, d$gen_idx, d$rep_idx,
    design_type = 1L
  )
  expect_equal(res$GMS, diag(res$MSG))
  expect_equal(res$EMS, diag(res$MSE))
})

# ==============================================================================
# NEW COVERAGE TESTS — targeting previously uncovered lines
# ==============================================================================

test_that("design_stats coerces non-numeric trait1 and trait2 (lines 100-101)", {
  d <- rcbd_data()

  # Pass as characters to trigger !is.numeric branch
  trait1_char <- as.character(d$trait1)
  trait2_char <- as.character(d$trait2)

  res <- selection.index:::design_stats(
    trait1 = trait1_char,
    trait2 = trait2_char,
    genotypes = d$gen_idx,
    replications = d$rep_idx,
    design_type = "RCBD",
    calc_type = "all"
  )

  expect_true(is.numeric(res$TSP))
  expect_true(is.finite(res$TSP))
})
