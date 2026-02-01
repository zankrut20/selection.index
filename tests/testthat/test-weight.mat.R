test_that("weight.mat converts dataframe to matrix correctly", {
  result <- weight.mat(data = weight)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), ncol(weight) - 1)
})

test_that("weight.mat removes first column", {
  result <- weight.mat(data = weight)
  
  # First column should be removed
  expect_false(colnames(weight)[1] %in% colnames(result))
})

test_that("weight.mat preserves column names", {
  result <- weight.mat(data = weight)
  
  # Check that remaining columns have correct names
  expected_names <- colnames(weight)[-1]
  expect_equal(colnames(result), expected_names)
})

test_that("weight.mat handles numeric conversion", {
  result <- weight.mat(data = weight)
  
  # All values should be numeric
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})

test_that("weight.mat preserves data values", {
  result <- weight.mat(data = weight)
  
  # Compare values (excluding first column)
  expected_matrix <- as.matrix(weight[, -1, drop = FALSE])
  expect_equal(result, expected_matrix)
})

test_that("weight.mat works with minimal dataframe", {
  # Create minimal test dataframe
  test_df <- data.frame(trait = c("T1", "T2"), w1 = c(1.0, 2.0))
  result <- weight.mat(data = test_df)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 1)
  expect_equal(nrow(result), 2)
})

test_that("weight.mat maintains row count", {
  result <- weight.mat(data = weight)
  
  expect_equal(nrow(result), nrow(weight))
})
