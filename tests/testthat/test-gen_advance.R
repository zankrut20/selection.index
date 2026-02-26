test_that("gen.advance calculates expected value", {
  seldata <- rnorm(2, 3, 4)
  gmat <- matrix(1:16, nrow = 4)
  pmat <- matrix(1:16, nrow = 4)
  weight <- matrix(1.00)
  GA <- round(gen_advance(phen_mat = pmat[1, 1], gen_mat = gmat[1, 1], weight_mat = weight), 2)
  result <- structure(2.060, .Dim = c(1L, 1L))
  expect_equal(result, GA)
})

test_that("gen.advance returns matrix output", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  GA <- gen_advance(phen_mat = pmat[1, 1], gen_mat = gmat[1, 1], weight_mat = weight[1, 2])

  expect_true(is.matrix(GA))
  expect_equal(dim(GA), c(1, 1))
})

test_that("gen.advance works with variance-covariance matrices", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

  # Test with single element
  GA1 <- gen_advance(phen_mat = pmat[1, 1], gen_mat = gmat[1, 1], weight_mat = weight[1, 2])
  expect_true(is.finite(GA1[1, 1]))
  expect_true(GA1[1, 1] > 0)
})

test_that("gen.advance handles matrix inputs correctly", {
  # Create symmetric positive definite matrices
  pmat <- matrix(c(4, 2, 2, 3), nrow = 2)
  gmat <- matrix(c(2, 1, 1, 1.5), nrow = 2)
  wmat <- matrix(c(1, 1), nrow = 2)

  GA <- gen_advance(phen_mat = pmat, gen_mat = gmat, weight_mat = wmat)
  expect_true(is.matrix(GA))
  expect_true(is.finite(GA[1, 1]))
})

test_that("gen.advance returns rounded values to 4 decimal places", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  GA <- gen_advance(phen_mat = pmat[1, 1], gen_mat = gmat[1, 1], weight_mat = weight[1, 2])

  # Check that value is rounded to 4 decimal places
  GA_str <- format(GA[1, 1], nsmall = 4)
  expect_true(nchar(sub(".*\\.", "", GA_str)) <= 4)
})
