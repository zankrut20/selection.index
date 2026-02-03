test_that("design.stats calculates sums of products correctly", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  trait2 <- seldata$dtf
  
  result <- design.stats(trait1, trait2, gen_idx, rep_idx, calc_type = "sums_of_products")
  
  # Check all expected components are present
  expect_true("CF" %in% names(result))
  expect_true("TSP" %in% names(result))
  expect_true("GSP" %in% names(result))
  expect_true("RSP" %in% names(result))
  expect_true("ESP" %in% names(result))
  
  # Check all values are numeric and finite
  expect_true(is.numeric(result$CF))
  expect_true(is.finite(result$CF))
  expect_true(all(is.finite(c(result$TSP, result$GSP, result$RSP, result$ESP))))
  
  # Check relationship: TSP = GSP + RSP + ESP
  expect_equal(result$TSP, result$GSP + result$RSP + result$ESP, tolerance = 1e-10)
})

test_that("design.stats calculates mean products correctly", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  
  result <- design.stats(trait1, trait1, gen_idx, rep_idx, calc_type = "mean_products")
  
  # Check expected components
  expect_true("GMP" %in% names(result))
  expect_true("EMP" %in% names(result))
  expect_true("DFG" %in% names(result))
  expect_true("DFE" %in% names(result))
  
  # Check all values are finite
  expect_true(all(is.finite(c(result$GMP, result$EMP))))
  
  # Check degrees of freedom are correct
  n_gen <- length(unique(gen_idx))
  n_rep <- length(unique(rep_idx))
  expect_equal(result$DFG, n_gen - 1)
  expect_equal(result$DFE, (n_gen - 1) * (n_rep - 1))
})

test_that("design.stats returns all components with calc_type='all'", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  trait2 <- seldata$dtf
  
  result <- design.stats(trait1, trait2, gen_idx, rep_idx, calc_type = "all")
  
  # Check all components are present
  expected_names <- c("CF", "TSP", "GSP", "RSP", "ESP", "GMP", "EMP", 
                      "DFG", "DFR", "DFE", "n_genotypes", "n_replications")
  expect_true(all(expected_names %in% names(result)))
  
  # Check all numeric values are finite
  numeric_vals <- c(result$CF, result$TSP, result$GSP, result$RSP, result$ESP,
                    result$GMP, result$EMP)
  expect_true(all(is.finite(numeric_vals)))
})

test_that("design.stats handles same trait (variance) correctly", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  
  # When trait1 = trait2, we're calculating variance components
  result <- design.stats(trait1, trait1, gen_idx, rep_idx, calc_type = "all")
  
  # For variance (same trait), TSP should be sum of squared deviations
  expect_true(result$TSP > 0)
  expect_true(result$GSP > 0)
  expect_true(result$ESP >= 0)  # Can be zero in extreme cases
  
  # GMP should be greater than or equal to EMP (genotype effect)
  expect_true(result$GMP >= result$EMP || abs(result$GMP - result$EMP) < 1e-10)
})

test_that("design.stats calculates anova_stats correctly", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  
  result <- design.stats(trait1, trait1, gen_idx, rep_idx, calc_type = "anova_stats")
  
  # Should only return degrees of freedom and basic stats
  expect_true("DFG" %in% names(result))
  expect_true("DFR" %in% names(result))
  expect_true("DFE" %in% names(result))
  expect_true("n_genotypes" %in% names(result))
  expect_true("n_replications" %in% names(result))
  
  # Should not include sums of products or mean products
  expect_false("TSP" %in% names(result))
  expect_false("GMP" %in% names(result))
})

test_that("design.stats degrees of freedom are correct", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  
  result <- design.stats(trait1, trait1, gen_idx, rep_idx, calc_type = "anova_stats")
  
  n_gen <- length(unique(gen_idx))
  n_rep <- length(unique(rep_idx))
  
  expect_equal(result$DFG, n_gen - 1)
  expect_equal(result$DFR, n_rep - 1)
  expect_equal(result$DFE, (n_gen - 1) * (n_rep - 1))
  expect_equal(result$n_genotypes, n_gen)
  expect_equal(result$n_replications, n_rep)
})

test_that("design.stats matches manual calculation", {
  # Create simple test data for manual verification
  gen_idx <- rep(1:3, each = 4)
  rep_idx <- rep(1:4, 3)
  trait1 <- c(10, 12, 11, 13,  # Gen 1
              15, 16, 14, 17,  # Gen 2
              20, 19, 21, 22)  # Gen 3
  trait2 <- c(5, 6, 7, 8,
              9, 10, 11, 12,
              13, 14, 15, 16)
  
  result <- design.stats(trait1, trait2, gen_idx, rep_idx, calc_type = "all")
  
  # Manual calculation of CF
  GT1 <- sum(trait1)
  GT2 <- sum(trait2)
  n_obs <- 3 * 4
  CF_manual <- (GT1 * GT2) / n_obs
  
  expect_equal(result$CF, CF_manual, tolerance = 1e-10)
  
  # Check TSP = sum(trait1 * trait2) - CF
  TSP_manual <- sum(trait1 * trait2) - CF_manual
  expect_equal(result$TSP, TSP_manual, tolerance = 1e-10)
})

test_that("design.stats works with different trait combinations", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  
  # Test with multiple trait pairs
  traits <- c("sypp", "dtf", "rpp")
  
  for(t1 in traits) {
    for(t2 in traits) {
      trait1 <- seldata[[t1]]
      trait2 <- seldata[[t2]]
      
      result <- design.stats(trait1, trait2, gen_idx, rep_idx, calc_type = "all")
      
      expect_true(all(is.finite(c(result$CF, result$TSP, result$GSP, 
                                  result$RSP, result$ESP, result$GMP, result$EMP))),
                  info = paste("Traits:", t1, "and", t2))
    }
  }
})

test_that("design.stats covariance matches gen.varcov", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  trait2 <- seldata$dtf
  
  # Use design engine
  design_result <- design.stats(trait1, trait2, gen_idx, rep_idx, calc_type = "mean_products")
  genetic_cov_engine <- (design_result$GMP - design_result$EMP) / design_result$n_replications
  
  # Use gen.varcov
  gmat <- gen.varcov(seldata[, c("sypp", "dtf")], 
                     seldata$treat, seldata$rep)
  genetic_cov_direct <- gmat[1, 2]
  
  # Should match within numerical tolerance
  expect_equal(genetic_cov_engine, genetic_cov_direct, tolerance = 1e-10)
})

test_that("design.stats default calc_type is 'all'", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  
  # Call without calc_type argument
  result <- design.stats(trait1, trait1, gen_idx, rep_idx)
  
  # Should have all components
  expected_names <- c("CF", "TSP", "GSP", "RSP", "ESP", "GMP", "EMP", 
                      "DFG", "DFR", "DFE", "n_genotypes", "n_replications")
  expect_true(all(expected_names %in% names(result)))
})

test_that("design.stats trait2 defaults to trait1", {
  gen_idx <- as.integer(as.factor(seldata$treat))
  rep_idx <- as.integer(as.factor(seldata$rep))
  trait1 <- seldata$sypp
  
  # Call without trait2 argument (should use trait1 for both)
  result1 <- design.stats(trait1, trait1, gen_idx, rep_idx, calc_type = "all")
  result2 <- design.stats(trait1, genotypes = gen_idx, replications = rep_idx, calc_type = "all")
  
  # Should produce identical results
  expect_equal(result1$CF, result2$CF)
  expect_equal(result1$TSP, result2$TSP)
  expect_equal(result1$GSP, result2$GSP)
})
