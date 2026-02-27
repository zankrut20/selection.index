test_that("gen.varcov and phen.varcov return correct dimensions", {
  gen <- gen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)
  phen <- phen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)

  expect_equal(nrow(gen), 7)
  expect_equal(ncol(gen), 7)
  expect_equal(nrow(phen), 7)
  expect_equal(ncol(phen), 7)
})

test_that("gen.varcov and phen.varcov return symmetric matrices", {
  gen <- gen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)
  phen <- phen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)

  expect_true(isSymmetric(gen))
  expect_true(isSymmetric(phen))
})

test_that("gen.varcov and phen.varcov have appropriate column names", {
  gen <- gen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)
  phen <- phen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)

  expect_equal(colnames(gen), colnames(seldata[, 3:9]))
  expect_equal(colnames(phen), colnames(seldata[, 3:9]))
  expect_equal(rownames(gen), colnames(seldata[, 3:9]))
  expect_equal(rownames(phen), colnames(seldata[, 3:9]))
})

test_that("gen.varcov returns values less than or equal to phen.varcov", {
  gen <- gen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)
  phen <- phen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)

  # Diagonal elements: genotypic variance should be <= phenotypic variance
  expect_true(all(diag(gen) <= diag(phen)))
})

test_that("gen.varcov and phen.varcov handle single trait", {
  gen_single <- gen_varcov(seldata[, 3, drop = FALSE], seldata$treat, seldata$rep)
  phen_single <- phen_varcov(seldata[, 3, drop = FALSE], seldata$treat, seldata$rep)

  expect_equal(dim(gen_single), c(1, 1))
  expect_equal(dim(phen_single), c(1, 1))
  expect_true(is.finite(gen_single[1, 1]))
  expect_true(is.finite(phen_single[1, 1]))
})

test_that("gen.varcov and phen.varcov work with different missing value methods", {
  # Create test data with missing values
  test_data <- seldata[, 3:9]
  test_data[1, 1] <- NA
  test_data[5, 3] <- NA

  # Test all methods
  methods <- c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")

  for (method in methods) {
    gen <- gen_varcov(test_data, seldata$treat, seldata$rep, method = method)
    phen <- phen_varcov(test_data, seldata$treat, seldata$rep, method = method)

    expect_equal(nrow(gen), 7, info = paste("Method:", method))
    expect_equal(nrow(phen), 7, info = paste("Method:", method))
    expect_true(all(is.finite(gen)), info = paste("Method:", method))
    expect_true(all(is.finite(phen)), info = paste("Method:", method))
  }
})

test_that("gen.varcov and phen.varcov return all finite values", {
  gen <- gen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)
  phen <- phen_varcov(seldata[, 3:9], seldata$treat, seldata$rep)

  expect_true(all(is.finite(gen)))
  expect_true(all(is.finite(phen)))
})

# ==============================================================================
# NEW COVERAGE TESTS — targeting previously uncovered lines in R/varcov.R
# ==============================================================================

# Helper to build minimal LSD and SPD datasets
.make_lsd <- function() {
  set.seed(123)
  data.frame(
    row_block = rep(1:3, each = 3),
    col_block = rep(1:3, 3),
    geno      = c(1, 2, 3, 2, 3, 1, 3, 1, 2),
    trait1    = rnorm(9, 30, 5),
    trait2    = rnorm(9, 20, 3)
  )
}

.make_spd <- function() {
  set.seed(123)
  data.frame(
    rep       = rep(1:3, each = 4),
    main_plot = rep(rep(1:2, each = 2), 3),
    sub_plot  = rep(1:2, 6),
    trait1    = rnorm(12, 50, 8),
    trait2    = rnorm(12, 10, 2)
  )
}

test_that(".calculate_varcov internal validations (lines 31-40)", {
  mat <- as.matrix(seldata[1:10, 3:4])
  gen_idx <- as.integer(as.factor(seldata$treat[1:10]))
  rep_idx <- as.integer(as.factor(seldata$rep[1:10]))

  # line 31: invalid design_type
  expect_error(
    selection.index:::.calculate_varcov(mat, gen_idx, rep_idx, design_type = 4L),
    "design_type must be 1 \\(RCBD\\), 2 \\(LSD\\), or 3 \\(SPD\\)."
  )

  # line 34: invalid cov_type
  expect_error(
    selection.index:::.calculate_varcov(mat, gen_idx, rep_idx, design_type = 1L, cov_type = 3L),
    "cov_type must be 1 \\(genotypic\\) or 2 \\(phenotypic\\)."
  )

  # line 37: LSD without col_idx
  expect_error(
    selection.index:::.calculate_varcov(mat, gen_idx, rep_idx, design_type = 2L, col_idx = NULL),
    "col_idx is required for Latin Square Design"
  )

  # line 40: SPD without main_idx
  expect_error(
    selection.index:::.calculate_varcov(mat, gen_idx, rep_idx, design_type = 3L, main_idx = NULL),
    "main_idx is required for Split Plot Design"
  )
})

test_that("gen_varcov validations for LSD and SPD (lines 126-137) and SPD formula (lines 67-70)", {
  lsd <- .make_lsd()
  spd <- .make_spd()

  # line 126: LSD without columns
  expect_error(
    gen_varcov(data = lsd[, c("trait1", "trait2")], genotypes = lsd$geno, replication = lsd$row_block, design_type = "LSD"),
    "Latin Square Design requires 'columns' parameter"
  )

  # line 129: LSD with columns
  res_lsd <- gen_varcov(data = lsd[, c("trait1", "trait2")], genotypes = lsd$geno, replication = lsd$row_block, columns = lsd$col_block, design_type = "LSD")
  expect_true(isSymmetric(res_lsd))

  # line 134: SPD without main_plots
  expect_error(
    gen_varcov(data = spd[, c("trait1", "trait2")], genotypes = spd$sub_plot, replication = spd$rep, design_type = "SPD"),
    "Split Plot Design requires 'main_plots' parameter"
  )

  # lines 137, 67-70: SPD with main_plots
  res_spd <- gen_varcov(data = spd[, c("trait1", "trait2")], genotypes = spd$sub_plot, replication = spd$rep, main_plots = spd$main_plot, design_type = "SPD")
  expect_true(isSymmetric(res_spd))
})

test_that("phen_varcov validations for LSD and SPD (lines 238-249)", {
  lsd <- .make_lsd()
  spd <- .make_spd()

  # line 238: LSD without columns
  expect_error(
    phen_varcov(data = lsd[, c("trait1", "trait2")], genotypes = lsd$geno, replication = lsd$row_block, design_type = "LSD"),
    "Latin Square Design requires 'columns' parameter"
  )

  # line 241: LSD with columns
  res_lsd <- phen_varcov(data = lsd[, c("trait1", "trait2")], genotypes = lsd$geno, replication = lsd$row_block, columns = lsd$col_block, design_type = "LSD")
  expect_true(isSymmetric(res_lsd))

  # line 246: SPD without main_plots
  expect_error(
    phen_varcov(data = spd[, c("trait1", "trait2")], genotypes = spd$sub_plot, replication = spd$rep, design_type = "SPD"),
    "Split Plot Design requires 'main_plots' parameter"
  )

  # line 249: SPD with main_plots
  res_spd <- phen_varcov(data = spd[, c("trait1", "trait2")], genotypes = spd$sub_plot, replication = spd$rep, main_plots = spd$main_plot, design_type = "SPD")
  expect_true(isSymmetric(res_spd))
})
