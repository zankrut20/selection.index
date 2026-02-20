# Comprehensive tests for C++ math primitives - Direct C++ function testing
# Goal: 100% coverage of src/math_primitives.cpp

# ==============================================================================
# TEST: cpp_grouped_sums (core function)
# ==============================================================================

test_that("cpp_grouped_sums handles basic grouping", {
  set.seed(123)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  group_idx <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L)
  
  result <- selection.index:::cpp_grouped_sums(data_mat, group_idx)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 3)
  
  # Verify sum for first group
  expect_equal(result[1, ], colSums(data_mat[1:3, ]))
})

test_that("cpp_grouped_sums handles empty data", {
  data_mat <- matrix(numeric(0), nrow = 0, ncol = 3)
  group_idx <- integer(0)
  
  result <- selection.index:::cpp_grouped_sums(data_mat, group_idx)
  expect_equal(dim(result), c(0, 3))
})

test_that("cpp_grouped_sums handles single group", {
  data_mat <- matrix(rnorm(15), nrow = 5, ncol = 3)
  group_idx <- rep(1L, 5)
  
  result <- selection.index:::cpp_grouped_sums(data_mat, group_idx)
  expect_equal(nrow(result), 1)
  expect_equal(result[1, ], colSums(data_mat))
})

test_that("cpp_grouped_sums validates group_idx size", {
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  group_idx <- c(1L, 2L, 3L)  # Wrong size
  
  expect_error(selection.index:::cpp_grouped_sums(data_mat, group_idx),
               "group_idx size must match")
})

test_that("cpp_grouped_sums validates positive integers", {
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  group_idx <- c(0L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L)  # Contains 0
  
  expect_error(selection.index:::cpp_grouped_sums(data_mat, group_idx),
               "must contain positive integers")
})

# ==============================================================================
# TEST: cpp_multi_grouped_sums
# ==============================================================================

test_that("cpp_multi_grouped_sums handles multiple groupings", {
  set.seed(456)
  data_mat <- matrix(rnorm(24), nrow = 8, ncol = 3)
  group_idx1 <- c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L)
  group_idx2 <- c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L)
  
  result <- selection.index:::cpp_multi_grouped_sums(data_mat, list(group_idx1, group_idx2))
  
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_equal(nrow(result[[1]]), 2)
  expect_equal(nrow(result[[2]]), 2)
})

test_that("cpp_multi_grouped_sums handles empty data", {
  data_mat <- matrix(numeric(0), nrow = 0, ncol = 3)
  group_indices <- list(integer(0), integer(0))
  
  result <- selection.index:::cpp_multi_grouped_sums(data_mat, group_indices)
  expect_equal(length(result), 2)
  expect_equal(dim(result[[1]]), c(0, 3))
})

test_that("cpp_multi_grouped_sums validates group sizes", {
  data_mat <- matrix(rnorm(24), nrow = 8, ncol = 3)
  group_idx1 <- c(1L, 1L, 2L, 2L)  # Wrong size
  group_idx2 <- c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L)
  
  expect_error(selection.index:::cpp_multi_grouped_sums(data_mat, list(group_idx1, group_idx2)),
               "size must match")
})

test_that("cpp_multi_grouped_sums validates positive integers", {
  data_mat <- matrix(rnorm(24), nrow = 8, ncol = 3)
  group_idx1 <- c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L)
  group_idx2 <- c(0L, 1L, 1L, 2L, 1L, 2L, 1L, 2L)  # Contains 0
  
  expect_error(selection.index:::cpp_multi_grouped_sums(data_mat, list(group_idx1, group_idx2)),
               "must contain positive integers")
})

# ==============================================================================
# TEST: cpp_crossprod_divided
# ==============================================================================

test_that("cpp_crossprod_divided computes correctly", {
  set.seed(789)
  sums1 <- matrix(rnorm(15), nrow = 5, ncol = 3)
  sums2 <- matrix(rnorm(15), nrow = 5, ncol = 3)
  divisor <- 10.0
  
  result <- selection.index:::cpp_crossprod_divided(sums1, sums2, divisor)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  
  # Verify computation: (t(sums1) %*% sums2) / divisor
  expected <- (t(sums1) %*% sums2) / divisor
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cpp_crossprod_divided handles identity case", {
  sums <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  result <- selection.index:::cpp_crossprod_divided(sums, sums, 1.0)
  
  expect_equal(result, diag(2), tolerance = 1e-10)
})

# ==============================================================================
# TEST: cpp_correction_factor_matrix
# ==============================================================================

test_that("cpp_correction_factor_matrix computes correctly", {
  set.seed(101)
  data_mat <- matrix(rnorm(40), nrow = 10, ncol = 4)
  
  result <- selection.index:::cpp_correction_factor_matrix(data_mat)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 4))
  expect_true(isSymmetric(result))
  
  # Verify computation
  grand_totals <- colSums(data_mat)
  expected <- outer(grand_totals, grand_totals) / nrow(data_mat)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cpp_correction_factor_matrix handles single column", {
  data_mat <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  result <- selection.index:::cpp_correction_factor_matrix(data_mat)
  
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1, 1], (15 * 15) / 5)
})

# ==============================================================================
# TEST: cpp_grand_means
# ==============================================================================

test_that("cpp_grand_means computes correctly", {
  set.seed(202)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  
  result <- selection.index:::cpp_grand_means(data_mat)
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 3)
  
  # Verify computation
  expected <- colMeans(data_mat)
  expect_equal(as.vector(result), expected, tolerance = 1e-10)
})

test_that("cpp_grand_means handles single row", {
  data_mat <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  result <- selection.index:::cpp_grand_means(data_mat)
  
  expect_equal(as.vector(result), c(1, 2, 3))
})

test_that("cpp_grand_means handles single column", {
  data_mat <- matrix(c(1, 2, 3, 4, 5), nrow = 5, ncol = 1)
  result <- selection.index:::cpp_grand_means(data_mat)
  
  expect_equal(as.vector(result), 3)
})

# ==============================================================================
# TEST: cpp_trait_minmax
# ==============================================================================

test_that("cpp_trait_minmax computes correctly", {
  set.seed(303)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  
  result <- selection.index:::cpp_trait_minmax(data_mat)
  
  expect_true(is.list(result))
  expect_true("min" %in% names(result))
  expect_true("max" %in% names(result))
  expect_equal(length(result$min), 3)
  expect_equal(length(result$max), 3)
  
  # Verify computation
  for (i in 1:3) {
    expect_equal(result$min[i], min(data_mat[, i]), tolerance = 1e-10)
    expect_equal(result$max[i], max(data_mat[, i]), tolerance = 1e-10)
  }
})

test_that("cpp_trait_minmax handles constant columns", {
  data_mat <- matrix(c(5, 5, 5, 1, 2, 3), nrow = 3, ncol = 2)
  result <- selection.index:::cpp_trait_minmax(data_mat)
  
  expect_equal(result$min[1], 5)
  expect_equal(result$max[1], 5)
  expect_equal(result$min[2], 1)
  expect_equal(result$max[2], 3)
})

# ==============================================================================
# TEST: cpp_genotype_means
# ==============================================================================

test_that("cpp_genotype_means computes correctly", {
  set.seed(404)
  data_mat <- matrix(rnorm(30), nrow = 12, ncol = 3)
  gen_idx <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L)
  
  result <- selection.index:::cpp_genotype_means(data_mat, gen_idx)
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 3)
  
  # Verify computation for first genotype
  expected_mean <- colMeans(data_mat[1:3, ])
  expect_equal(result[1, ], expected_mean, tolerance = 1e-10)
})

test_that("cpp_genotype_means handles unbalanced groups", {
  data_mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 9, ncol = 1)
  gen_idx <- c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L)
  
  result <- selection.index:::cpp_genotype_means(data_mat, gen_idx)
  
  expect_equal(result[1, 1], 1.5)  # mean(1, 2)
  expect_equal(result[2, 1], 4.0)  # mean(3, 4, 5)
  expect_equal(result[3, 1], 7.5)  # mean(6, 7, 8, 9)
})

# ==============================================================================
# TEST: cpp_extract_submatrix
# ==============================================================================

test_that("cpp_extract_submatrix extracts correctly", {
  mat <- matrix(1:16, nrow = 4, ncol = 4)
  indices <- c(1L, 3L)  # Extract rows/columns 1 and 3
  
  result <- selection.index:::cpp_extract_submatrix(mat, indices)
  
  expect_equal(dim(result), c(2, 2))
  expect_equal(result[1, 1], mat[1, 1])
  expect_equal(result[1, 2], mat[1, 3])
  expect_equal(result[2, 1], mat[3, 1])
  expect_equal(result[2, 2], mat[3, 3])
})

test_that("cpp_extract_submatrix handles single element", {
  mat <- matrix(1:16, nrow = 4, ncol = 4)
  indices <- c(2L)
  
  result <- selection.index:::cpp_extract_submatrix(mat, indices)
  
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1, 1], mat[2, 2])
})

test_that("cpp_extract_submatrix handles full matrix", {
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  indices <- c(1L, 2L, 3L)
  
  result <- selection.index:::cpp_extract_submatrix(mat, indices)
  
  expect_equal(result, mat)
})

test_that("cpp_extract_submatrix handles reordering", {
  mat <- matrix(1:9, nrow = 3, ncol = 3)
  indices <- c(3L, 1L)  # Reverse order
  
  result <- selection.index:::cpp_extract_submatrix(mat, indices)
  
  expect_equal(result[1, 1], mat[3, 3])
  expect_equal(result[1, 2], mat[3, 1])
  expect_equal(result[2, 1], mat[1, 3])
  expect_equal(result[2, 2], mat[1, 1])
})

# ==============================================================================
# TEST: cpp_extract_vector
# ==============================================================================

test_that("cpp_extract_vector extracts correctly", {
  mat <- matrix(1:12, nrow = 4, ncol = 3)
  row_indices <- c(1L, 3L)
  col_index <- 1L  # Second column (0-based in C++)
  
  result <- selection.index:::cpp_extract_vector(mat, row_indices, col_index)
  
  expect_equal(length(result), 2)
  expect_equal(result[1], mat[1, 2])
  expect_equal(result[2], mat[3, 2])
})

test_that("cpp_extract_vector handles single element", {
  mat <- matrix(1:12, nrow = 4, ncol = 3)
  row_indices <- c(2L)
  col_index <- 2L
  
  result <- selection.index:::cpp_extract_vector(mat, row_indices, col_index)
  
  expect_equal(length(result), 1)
  expect_equal(result[1], mat[2, 3])
})

test_that("cpp_extract_vector handles all rows", {
  mat <- matrix(1:12, nrow = 4, ncol = 3)
  row_indices <- c(1L, 2L, 3L, 4L)
  col_index <- 0L
  
  result <- selection.index:::cpp_extract_vector(mat, row_indices, col_index)
  
  expect_equal(as.vector(result), mat[, 1])
})

# ==============================================================================
# TEST: cpp_symmetric_solve
# ==============================================================================

test_that("cpp_symmetric_solve solves linear system", {
  set.seed(505)
  A <- matrix(c(4, 1, 1, 3), nrow = 2, ncol = 2)
  b <- c(1, 2)
  
  result <- selection.index:::cpp_symmetric_solve(A, b)
  
  expect_equal(length(result), 2)
  
  # Verify solution: A %*% result should equal b
  verification <- A %*% result
  expect_equal(as.vector(verification), b, tolerance = 1e-10)
})

test_that("cpp_symmetric_solve handles identity matrix", {
  A <- diag(3)
  b <- c(1, 2, 3)
  
  result <- selection.index:::cpp_symmetric_solve(A, b)
  
  expect_equal(as.vector(result), b, tolerance = 1e-10)
})

test_that("cpp_symmetric_solve handles larger systems", {
  set.seed(606)
  # Create symmetric positive definite matrix
  M <- matrix(rnorm(16), 4, 4)
  A <- t(M) %*% M  # Guaranteed positive definite
  b <- rnorm(4)
  
  result <- selection.index:::cpp_symmetric_solve(A, b)
  
  # Verify solution
  verification <- A %*% result
  expect_equal(as.vector(verification), b, tolerance = 1e-9)
})

# ==============================================================================
# TEST: cpp_quadratic_form
# ==============================================================================

test_that("cpp_quadratic_form computes correctly", {
  x <- c(1, 2)
  A <- matrix(c(3, 1, 1, 2), nrow = 2, ncol = 2)
  y <- c(2, 1)
  
  result <- selection.index:::cpp_quadratic_form(x, A, y)
  
  # Verify: t(x) %*% A %*% y
  expected <- as.numeric(t(x) %*% A %*% y)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cpp_quadratic_form handles identity matrix", {
  x <- c(1, 2, 3)
  A <- diag(3)
  y <- c(4, 5, 6)
  
  result <- selection.index:::cpp_quadratic_form(x, A, y)
  
  # With identity matrix: result should be dot product
  expected <- sum(x * y)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cpp_quadratic_form handles zero vectors", {
  x <- c(0, 0)
  A <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  y <- c(5, 6)
  
  result <- selection.index:::cpp_quadratic_form(x, A, y)
  
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("cpp_quadratic_form handles larger matrices", {
  set.seed(707)
  x <- rnorm(5)
  A <- matrix(rnorm(25), nrow = 5, ncol = 5)
  y <- rnorm(5)
  
  result <- selection.index:::cpp_quadratic_form(x, A, y)
  expected <- as.numeric(t(x) %*% A %*% y)
  
  expect_equal(result, expected, tolerance = 1e-9)
})

# ==============================================================================
# TEST: cpp_quadratic_form_sym
# ==============================================================================

test_that("cpp_quadratic_form_sym computes correctly", {
  x <- c(1, 2)
  A <- matrix(c(3, 1, 1, 2), nrow = 2, ncol = 2)
  
  result <- selection.index:::cpp_quadratic_form_sym(x, A)
  
  # Verify: t(x) %*% A %*% x
  expected <- as.numeric(t(x) %*% A %*% x)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cpp_quadratic_form_sym handles identity matrix", {
  x <- c(1, 2, 3)
  A <- diag(3)
  
  result <- selection.index:::cpp_quadratic_form_sym(x, A)
  
  # With identity matrix: result should be sum of squares
  expected <- sum(x^2)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("cpp_quadratic_form_sym handles zero vector", {
  x <- c(0, 0, 0)
  A <- matrix(rnorm(9), nrow = 3, ncol = 3)
  
  result <- selection.index:::cpp_quadratic_form_sym(x, A)
  
  expect_equal(result, 0, tolerance = 1e-10)
})

test_that("cpp_quadratic_form_sym handles larger matrices", {
  set.seed(808)
  x <- rnorm(6)
  # Create symmetric matrix
  M <- matrix(rnorm(36), 6, 6)
  A <- (M + t(M)) / 2
  
  result <- selection.index:::cpp_quadratic_form_sym(x, A)
  expected <- as.numeric(t(x) %*% A %*% x)
  
  expect_equal(result, expected, tolerance = 1e-9)
})

# ==============================================================================
# TEST: cpp_correction_factor (alternative implementation)
# ==============================================================================

test_that("cpp_correction_factor computes correctly", {
  total_sums <- c(10, 20, 30)
  n_obs <- 50L
  
  result <- selection.index:::cpp_correction_factor(total_sums, n_obs)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 3))
  expect_true(isSymmetric(result))
  
  # Verify computation
  expect_equal(result[1, 1], (10 * 10) / 50)
  expect_equal(result[1, 2], (10 * 20) / 50)
  expect_equal(result[2, 3], (20 * 30) / 50)
})

test_that("cpp_correction_factor handles single trait", {
  total_sums <- c(100)
  n_obs <- 25L
  
  result <- selection.index:::cpp_correction_factor(total_sums, n_obs)
  
  expect_equal(dim(result), c(1, 1))
  expect_equal(result[1, 1], (100 * 100) / 25)
})

# ==============================================================================
# TEST: cpp_total_sum_of_products
# ==============================================================================

test_that("cpp_total_sum_of_products computes correctly", {
  set.seed(909)
  data_mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
  total_sums <- colSums(data_mat)
  CF <- selection.index:::cpp_correction_factor(total_sums, nrow(data_mat))
  
  result <- selection.index:::cpp_total_sum_of_products(data_mat, CF)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_true(isSymmetric(result))
  
  # Verify computation manually
  TSP_manual <- matrix(0, 2, 2)
  for (i in 1:2) {
    for (j in 1:2) {
      TSP_manual[i, j] <- sum(data_mat[, i] * data_mat[, j]) - CF[i, j]
    }
  }
  expect_equal(result, TSP_manual, tolerance = 1e-10)
})

test_that("cpp_total_sum_of_products matches R implementation", {
  set.seed(1010)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  total_sums <- colSums(data_mat)
  CF <- selection.index:::cpp_correction_factor(total_sums, nrow(data_mat))
  
  result_cpp <- selection.index:::cpp_total_sum_of_products(data_mat, CF)
  
  # R implementation
  result_r <- t(data_mat) %*% data_mat - CF
  
  expect_equal(result_cpp, result_r, tolerance = 1e-10)
})

# ==============================================================================
# TEST: cpp_grouped_sum_of_products
# ==============================================================================

test_that("cpp_grouped_sum_of_products computes correctly", {
  set.seed(1111)
  data_mat <- matrix(rnorm(24), nrow = 12, ncol = 2)
  group_idx <- rep(1:4, each = 3)
  
  group_sums <- selection.index:::cpp_grouped_sums(data_mat, group_idx)
  group_counts <- as.integer(table(group_idx))
  total_sums <- colSums(data_mat)
  CF <- selection.index:::cpp_correction_factor(total_sums, nrow(data_mat))
  
  result <- selection.index:::cpp_grouped_sum_of_products(group_sums, group_counts, CF)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_true(isSymmetric(result))
})

test_that("cpp_grouped_sum_of_products handles single group", {
  data_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  group_sums <- matrix(colSums(data_mat), nrow = 1)
  group_counts <- 3L
  CF <- matrix(c(36/3, 45/3, 45/3, 63/3), nrow = 2)
  
  result <- selection.index:::cpp_grouped_sum_of_products(group_sums, group_counts, CF)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
})

# ==============================================================================
# TEST: cpp_mean_squares
# ==============================================================================

test_that("cpp_mean_squares computes correctly", {
  SP <- matrix(c(100, 50, 50, 80), nrow = 2, ncol = 2)
  df <- 10L
  
  result <- selection.index:::cpp_mean_squares(SP, df)
  
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2, 2))
  expect_equal(result, SP / df, tolerance = 1e-10)
})

test_that("cpp_mean_squares handles single degree of freedom", {
  SP <- matrix(c(50, 25, 25, 40), nrow = 2, ncol = 2)
  df <- 1L
  
  result <- selection.index:::cpp_mean_squares(SP, df)
  
  expect_equal(result, SP, tolerance = 1e-10)
})

test_that("cpp_mean_squares handles larger matrices", {
  set.seed(1212)
  SP <- matrix(rnorm(16), nrow = 4, ncol = 4)
  SP <- (SP + t(SP)) / 2  # Make symmetric
  df <- 20L
  
  result <- selection.index:::cpp_mean_squares(SP, df)
  
  expect_equal(result, SP / df, tolerance = 1e-10)
  expect_true(isSymmetric(result))
})
