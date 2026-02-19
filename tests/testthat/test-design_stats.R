# Comprehensive test for fixed design_stats.R - testing actual functionality

test_that("design_stats RCBD with fixed C++ code", {
  # Use realistic test data that previously caused crashes
  trait1 <- c(10, 12, 11, 13, 14, 12, 15, 16, 14)
  trait2 <- c(15, 16, 14, 18, 17, 15, 19, 20, 18)
  gen_idx <- c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L)  # 3 genotypes
  rep_idx <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)  # 3 reps
  
  # Test all calc_types that previously crashed
  result1 <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx, 
                                           design_type = "RCBD", calc_type = "sums_of_products")
  expect_true("CF" %in% names(result1))
  expect_true("TSP" %in% names(result1))
  expect_true("GSP" %in% names(result1))
  expect_true("RSP" %in% names(result1))
  expect_true("ESP" %in% names(result1))
  expect_true(all(is.finite(c(result1$CF, result1$TSP, result1$GSP, result1$RSP, result1$ESP))))
  
  result2 <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx, 
                                           design_type = "RCBD", calc_type = "mean_products")
  expect_true("GMP" %in% names(result2))
  expect_true("EMP" %in% names(result2))
  expect_true(all(is.finite(c(result2$GMP, result2$EMP))))
  
  result3 <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx, 
                                           design_type = "RCBD", calc_type = "all")
  expect_true(all(c("CF", "TSP", "GSP", "RSP", "ESP", "GMP", "EMP", "DFG", "DFR", "DFE") %in% names(result3)))
})

test_that("design_stats LSD with fixed C++ code", {
  trait1 <- c(10, 12, 11, 13, 14, 12, 15, 16, 14)
  trait2 <- c(15, 16, 14, 18, 17, 15, 19, 20, 18)
  gen_idx <- c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L)
  rep_idx <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L)
  columns <- c(1L, 2L, 3L, 2L, 3L, 1L, 3L, 1L, 2L)  # Latin square pattern
  
  result <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx, 
                                         columns = columns, design_type = "LSD", 
                                         calc_type = "all")
  expect_true("CSP" %in% names(result))  # Column sum of products
  expect_equal(result$design_type, "LSD")
})

test_that("design_stats SPD with fixed C++ code", {
  trait1 <- c(10, 12, 11, 13, 14, 12)
  trait2 <- c(15, 16, 14, 18, 17, 15)
  gen_idx <- c(1L, 2L, 1L, 2L, 1L, 2L)  # Sub-plot treatments  
  rep_idx <- c(1L, 1L, 2L, 2L, 3L, 3L)  # Blocks
  main_plots <- c(1L, 1L, 1L, 1L, 2L, 2L)  # Main plot treatments
  
  # Ensure all vectors are same length and correct types
  expect_equal(length(trait1), length(trait2))
  expect_equal(length(trait1), length(gen_idx))
  expect_equal(length(trait1), length(rep_idx))
  expect_equal(length(trait1), length(main_plots))
  
  result <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx,
                                         main_plots = main_plots, design_type = "SPD",
                                         calc_type = "all")
  expect_true("MSP" %in% names(result))     # Main plot sum of products
  expect_true("ESP_MAIN" %in% names(result)) # Main plot error
  expect_equal(result$design_type, "SPD")
  expect_equal(result$n_main_plots, 2)   # Should detect 2 main plot treatments
})

test_that("design_stats_api with fixed C++ code", {
  trait1 <- c(10, 12, 11, 13)
  trait2 <- c(15, 16, 14, 18)
  data_mat <- cbind(trait1, trait2)
  gen_idx <- c(1L, 2L, 1L, 2L)
  rep_idx <- c(1L, 1L, 2L, 2L)
  
  # This should now work without crashes
  result <- selection.index:::design_stats_api(data_mat, gen_idx, rep_idx, design_type = 1L)
  expect_true("MSG" %in% names(result))
  expect_true("MSE" %in% names(result))
  expect_equal(dim(result$MSG), c(2, 2))
  expect_equal(dim(result$MSE), c(2, 2))
})

test_that("design_stats with seldata subset", {
  # Test with actual seldata that was causing crashes
  data(seldata)
  idx <- which(seldata[["treat"]] %in% c("G1", "G2", "G3"))
  trait1 <- as.numeric(seldata[["sypp"]][idx])
  trait2 <- as.numeric(seldata[["dtf"]][idx])
  gen_idx <- as.integer(as.factor(seldata[["treat"]][idx]))
  rep_idx <- as.integer(seldata[["rep"]][idx])
  
  # This should now work without matrix dimension errors  
  result <- selection.index:::design_stats(trait1, trait2, gen_idx, rep_idx,
                                          design_type = "RCBD", calc_type = "all")
  expect_equal(result$n_genotypes, 3)
  expect_equal(result$n_replications, 3)
  expect_true(all(is.finite(c(result$CF, result$GMP, result$EMP))))
})