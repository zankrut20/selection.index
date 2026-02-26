test_that("lpsi basic functionality works", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  result <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 7) # 7 traits, choose 1
  expect_true("ID" %in% colnames(result))
  expect_true("GA" %in% colnames(result))
  expect_true("PRE" %in% colnames(result))
  expect_true("Rank" %in% colnames(result))
})

test_that("lpsi excluding_trait with numeric indices works", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Exclude trait 1 (sypp)
  result <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = 1
  )

  # Should have choose(6, 2) = 15 combinations (excluding trait 1)
  expect_equal(nrow(result), 15)

  # No combination should contain trait 1
  ids <- strsplit(result$ID, ", ")
  has_trait_1 <- sapply(ids, function(x) "1" %in% x)
  expect_false(any(has_trait_1))
})

test_that("lpsi excluding_trait with multiple numeric indices works", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Exclude traits 1 and 2
  result <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = c(1, 2)
  )

  # Should have choose(5, 2) = 10 combinations
  expect_equal(nrow(result), 10)

  # No combination should contain trait 1 or 2
  ids <- strsplit(result$ID, ", ")
  has_excluded <- sapply(ids, function(x) any(c("1", "2") %in% x))
  expect_false(any(has_excluded))
})

test_that("lpsi excluding_trait with character trait names works", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Matrices should have column names
  expect_false(is.null(colnames(pmat)))
  expect_equal(colnames(pmat)[1], "sypp")

  # Exclude "sypp" by name
  result <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = "sypp"
  )

  # Should have choose(6, 2) = 15 combinations (excluding sypp)
  expect_equal(nrow(result), 15)

  # No combination should contain trait 1 (sypp)
  ids <- strsplit(result$ID, ", ")
  has_trait_1 <- sapply(ids, function(x) "1" %in% x)
  expect_false(any(has_trait_1))
})

test_that("lpsi excluding_trait with multiple character names works", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Exclude "sypp" and "dtf" by name
  result <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = c("sypp", "dtf")
  )

  # Should have choose(5, 2) = 10 combinations
  expect_equal(nrow(result), 10)

  # No combination should contain trait 1 or 2
  ids <- strsplit(result$ID, ", ")
  has_excluded <- sapply(ids, function(x) any(c("1", "2") %in% x))
  expect_false(any(has_excluded))
})

test_that("lpsi excluding_trait with data frame columns works", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Exclude using data columns (seldata[,3:4] = sypp, dtf)
  result <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = seldata[, 3:4]
  )

  # Should have choose(5, 2) = 10 combinations
  expect_equal(nrow(result), 10)

  # No combination should contain trait 1 or 2
  ids <- strsplit(result$ID, ", ")
  has_excluded <- sapply(ids, function(x) any(c("1", "2") %in% x))
  expect_false(any(has_excluded))
})

test_that("lpsi excluding all traits returns empty data frame", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Exclude all 7 traits
  result <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = 1:7
  )

  expect_equal(nrow(result), 0)
  expect_true(is.data.frame(result))
  expect_true("ID" %in% colnames(result))
})

test_that("lpsi with GAY calculates PRE correctly", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  result <- lpsi(
    ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, GAY = 1.075
  )

  expect_true("PRE" %in% colnames(result))
  expect_true(all(is.finite(result$PRE)))
})

test_that("lpsi warning for invalid character trait names", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  expect_warning(
    lpsi(
      ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
      wcol = 1, excluding_trait = "invalid_trait"
    ),
    "None of the specified trait names found"
  )
})

test_that("lpsi error for character names without pmat colnames", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Remove column names
  pmat_no_names <- pmat
  colnames(pmat_no_names) <- NULL

  expect_error(
    lpsi(
      ncomb = 2, pmat = pmat_no_names, gmat = gmat, wmat = weight[, -1],
      wcol = 1, excluding_trait = "sypp"
    ),
    "pmat must have column names"
  )
})

test_that("lpsi error for data frame without column names", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Create data frame without column names
  bad_df <- data.frame(matrix(1:10, ncol = 2))
  colnames(bad_df) <- NULL

  expect_error(
    lpsi(
      ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
      wcol = 1, excluding_trait = bad_df
    ),
    "excluding_trait data must have column names"
  )
})

test_that("lpsi excluding_trait produces correct combination count", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Without exclusion: choose(7, 2) = 21
  result_all <- lpsi(ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)
  expect_equal(nrow(result_all), 21)

  # Exclude 1 trait: choose(6, 2) = 15
  result_excl1 <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = 1
  )
  expect_equal(nrow(result_excl1), 15)

  # Exclude 2 traits: choose(5, 2) = 10
  result_excl2 <- lpsi(
    ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1],
    wcol = 1, excluding_trait = c(1, 2)
  )
  expect_equal(nrow(result_excl2), 10)
})

test_that("lpsi handles different wcol values", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  result1 <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)
  result2 <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 2)

  # Different weight columns should give different results
  expect_false(identical(result1$GA, result2$GA))
})

test_that("lpsi returns proper metrics", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  result <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  # Check all required metrics exist
  expect_true("Delta_G" %in% colnames(result))
  expect_true("rHI" %in% colnames(result))
  expect_true("hI2" %in% colnames(result))

  # Check metric ranges
  expect_true(all(result$hI2 >= 0 & result$hI2 <= 1))
  expect_true(all(result$rHI >= 0 & result$rHI <= 1))
})
