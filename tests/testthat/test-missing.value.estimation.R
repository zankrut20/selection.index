# Test the modular missing.value.estimation function

test_that("missing.value.estimation works with all three methods", {
  # Create sample RCBD data with missing values
  set.seed(123)
  n_gen <- 5
  n_rep <- 3
  n_traits <- 4
  n_obs <- n_gen * n_rep
  
  # Generate test data
  data_mat <- matrix(rnorm(n_obs * n_traits), nrow = n_obs, ncol = n_traits)
  colnames(data_mat) <- paste0("trait", 1:n_traits)
  
  gen_idx <- rep(1:n_gen, each = n_rep)
  rep_idx <- rep(1:n_rep, times = n_gen)
  
  # Introduce missing values
  missing_positions <- c(3, 7, 12, 15)
  data_with_missing <- data_mat
  data_with_missing[missing_positions, 1] <- NA
  data_with_missing[c(5, 10), 2] <- NA
  
  # Test REML method
  result_reml <- missing.value.estimation(data_with_missing, gen_idx, rep_idx, method = "REML")
  expect_true(all(is.finite(result_reml)))
  expect_equal(dim(result_reml), dim(data_with_missing))
  expect_false(any(is.na(result_reml)))
  
  # Test Yates method
  result_yates <- missing.value.estimation(data_with_missing, gen_idx, rep_idx, method = "Yates")
  expect_true(all(is.finite(result_yates)))
  expect_equal(dim(result_yates), dim(data_with_missing))
  expect_false(any(is.na(result_yates)))
  
  # Test Healy & Westmacott method
  result_healy <- missing.value.estimation(data_with_missing, gen_idx, rep_idx, method = "Healy")
  expect_true(all(is.finite(result_healy)))
  expect_equal(dim(result_healy), dim(data_with_missing))
  expect_false(any(is.na(result_healy)))
  
  # Test Regression method
  result_regression <- missing.value.estimation(data_with_missing, gen_idx, rep_idx, method = "Regression")
  expect_true(all(is.finite(result_regression)))
  expect_equal(dim(result_regression), dim(data_with_missing))
  expect_false(any(is.na(result_regression)))
  
  # Test Mean substitution method
  result_mean <- missing.value.estimation(data_with_missing, gen_idx, rep_idx, method = "Mean")
  expect_true(all(is.finite(result_mean)))
  expect_equal(dim(result_mean), dim(data_with_missing))
  expect_false(any(is.na(result_mean)))
  
  # Test Bartlett method
  result_bartlett <- missing.value.estimation(data_with_missing, gen_idx, rep_idx, method = "Bartlett")
  expect_true(all(is.finite(result_bartlett)))
  expect_equal(dim(result_bartlett), dim(data_with_missing))
  expect_false(any(is.na(result_bartlett)))
  
  # Test that non-missing values are preserved
  non_missing_idx <- which(is.finite(data_with_missing[, 3]))
  expect_equal(result_reml[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_yates[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_healy[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_regression[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_mean[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
  expect_equal(result_bartlett[non_missing_idx, 3], data_with_missing[non_missing_idx, 3])
})

test_that("missing.value.estimation returns data unchanged when no missing values", {
  set.seed(456)
  data_mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  gen_idx <- rep(1:5, each = 2)
  rep_idx <- rep(1:2, times = 5)
  
  result <- missing.value.estimation(data_mat, gen_idx, rep_idx)
  expect_equal(result, data_mat)
})

test_that("gen.varcov and phen.varcov use missing.value.estimation correctly", {
  # Load test data
  data(seldata, package = "selection.index")
  
  # Create data with missing values
  test_data <- seldata[, 3:9]
  test_data[3, 1] <- NA
  test_data[7, 2] <- NA
  test_data[15, 3] <- NA
  
  # Test that functions work with all three methods
  gen_reml <- gen.varcov(data = test_data, 
                         genotypes = seldata$treat, 
                         replication = seldata$rep,
                         method = "REML")
  expect_true(is.matrix(gen_reml))
  expect_equal(nrow(gen_reml), ncol(test_data))
  
  gen_yates <- gen.varcov(data = test_data, 
                          genotypes = seldata$treat, 
                          replication = seldata$rep,
                          method = "Yates")
  expect_true(is.matrix(gen_yates))
  expect_equal(nrow(gen_yates), ncol(test_data))
  
  gen_healy <- gen.varcov(data = test_data, 
                          genotypes = seldata$treat, 
                          replication = seldata$rep,
                          method = "Healy")
  expect_true(is.matrix(gen_healy))
  expect_equal(nrow(gen_healy), ncol(test_data))
  
  gen_regression <- gen.varcov(data = test_data, 
                               genotypes = seldata$treat, 
                               replication = seldata$rep,
                               method = "Regression")
  expect_true(is.matrix(gen_regression))
  expect_equal(nrow(gen_regression), ncol(test_data))
  
  gen_mean <- gen.varcov(data = test_data, 
                         genotypes = seldata$treat, 
                         replication = seldata$rep,
                         method = "Mean")
  expect_true(is.matrix(gen_mean))
  expect_equal(nrow(gen_mean), ncol(test_data))
  
  gen_bartlett <- gen.varcov(data = test_data, 
                             genotypes = seldata$treat, 
                             replication = seldata$rep,
                             method = "Bartlett")
  expect_true(is.matrix(gen_bartlett))
  expect_equal(nrow(gen_bartlett), ncol(test_data))
  
  phen_reml <- phen.varcov(data = test_data, 
                           genotypes = seldata$treat, 
                           replication = seldata$rep,
                           method = "REML")
  expect_true(is.matrix(phen_reml))
  expect_equal(nrow(phen_reml), ncol(test_data))
  
  phen_yates <- phen.varcov(data = test_data, 
                            genotypes = seldata$treat, 
                            replication = seldata$rep,
                            method = "Yates")
  expect_true(is.matrix(phen_yates))
  expect_equal(nrow(phen_yates), ncol(test_data))
  
  phen_healy <- phen.varcov(data = test_data, 
                            genotypes = seldata$treat, 
                            replication = seldata$rep,
                            method = "Healy")
  expect_true(is.matrix(phen_healy))
  expect_equal(nrow(phen_healy), ncol(test_data))
  
  phen_regression <- phen.varcov(data = test_data, 
                                 genotypes = seldata$treat, 
                                 replication = seldata$rep,
                                 method = "Regression")
  expect_true(is.matrix(phen_regression))
  expect_equal(nrow(phen_regression), ncol(test_data))
  
  phen_mean <- phen.varcov(data = test_data, 
                           genotypes = seldata$treat, 
                           replication = seldata$rep,
                           method = "Mean")
  expect_true(is.matrix(phen_mean))
  expect_equal(nrow(phen_mean), ncol(test_data))
  
  phen_bartlett <- phen.varcov(data = test_data, 
                               genotypes = seldata$treat, 
                               replication = seldata$rep,
                               method = "Bartlett")
  expect_true(is.matrix(phen_bartlett))
  expect_equal(nrow(phen_bartlett), ncol(test_data))
})

test_that("functions warn when missing values present but method not specified", {
  # Load test data
  data(seldata, package = "selection.index")
  
  # Create data with missing values
  test_data <- seldata[, 3:9]
  test_data[3, 1] <- NA
  test_data[7, 2] <- NA
  
  # Test gen.varcov warns when method not provided
  expect_warning(
    gen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep),
    "Missing values detected in data. Using default method 'REML' for imputation."
  )
  
  # Test phen.varcov warns when method not provided
  expect_warning(
    phen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep),
    "Missing values detected in data. Using default method 'REML' for imputation."
  )
  
  # Test no warning when method is explicitly provided
  expect_silent(
    gen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep, method = "REML")
  )
  
  expect_silent(
    phen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep, method = "Mean")
  )
})

test_that("functions work without warning when no missing values", {
  # Load test data
  data(seldata, package = "selection.index")
  
  # Use complete data (no missing values)
  test_data <- seldata[, 3:9]
  
  # Should not warn even without method specified
  expect_silent(
    gen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  )
  
  expect_silent(
    phen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  )
  
  # Verify matrices are returned
  result_gen <- gen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  result_phen <- phen.varcov(data = test_data, genotypes = seldata$treat, replication = seldata$rep)
  
  expect_true(is.matrix(result_gen))
  expect_true(is.matrix(result_phen))
  expect_equal(nrow(result_gen), ncol(test_data))
  expect_equal(nrow(result_phen), ncol(test_data))
})
