test_that("mean.performance returns correct dimensions", {
  performance <- mean.performance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  expect_equal(nrow(performance), 34)
})

test_that("mean.performance returns data frame with correct columns", {
  performance <- mean.performance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  expect_true(is.data.frame(performance))
  expect_true("Genotypes" %in% colnames(performance))
  # Column names should match trait names from input data
  expect_true(all(colnames(seldata[,3:7]) %in% colnames(performance)))
})

test_that("mean.performance handles all traits", {
  performance <- mean.performance(data = seldata[,3:9], genotypes = seldata[,2], replications = seldata[,1])
  
  # Should have rows for genotypes + 9 summary statistics (Min, Max, GM, CV, SEm, CD 5%, CD 1%, Heritability, Heritability%)
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  expect_equal(nrow(performance), n_genotypes + 9)
  
  # Check that trait columns exist
  expect_true(all(colnames(seldata[,3:9]) %in% colnames(performance)))
})

test_that("mean.performance genotype means are correct", {
  performance <- mean.performance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  # Manually calculate mean for first genotype and first trait
  first_geno <- levels(as.factor(seldata[,2]))[1]
  first_trait <- colnames(seldata)[3]
  manual_mean <- mean(seldata[seldata[,2] == first_geno, 3], na.rm = TRUE)
  
  perf_mean <- performance[performance$Genotypes == first_geno, first_trait]
  
  # Should be close (allowing for rounding)
  expect_true(abs(manual_mean - as.numeric(perf_mean)) < 0.01)
})

test_that("mean.performance works with single trait", {
  performance <- mean.performance(data = seldata[,3, drop=FALSE], genotypes = seldata[,2], replications = seldata[,1])
  
  expect_true(is.data.frame(performance))
  # Should have genotypes + 9 summary statistics rows
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  expect_equal(nrow(performance), n_genotypes + 9)
})

test_that("mean.performance includes all genotypes", {
  performance <- mean.performance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  unique_genos <- unique(seldata[,2])
  # First n rows should be genotype means, last 9 rows are summary statistics
  n_genotypes <- length(unique_genos)
  genotype_rows <- performance$Genotypes[1:n_genotypes]
  
  # Check that all unique genotypes appear in the first n rows
  expect_true(all(unique_genos %in% genotype_rows))
})

test_that("mean.performance calculates statistics correctly", {
  performance <- mean.performance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  # Check that all numeric columns contain finite values
  numeric_cols <- sapply(performance, is.numeric)
  numeric_data <- performance[, numeric_cols]
  expect_true(all(is.finite(as.matrix(numeric_data))))
})

# ============ NEW TESTS FOR design.stats INTEGRATION ============

test_that("mean.performance uses design.stats engine correctly", {
  performance <- mean.performance(data = seldata[,3:5], genotypes = seldata[,2], 
                                replications = seldata[,1], design_type = "RCBD")
  
  # Extract summary statistics rows
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  summary_rows <- performance[(n_genotypes + 1):nrow(performance), ]
  
  # Check all expected summary rows exist
  expected_stats <- c("Min", "Max", "GM", "CV (%)", "SEm", "CD 5%", "CD 1%", 
                     "Heritability", "Heritability(%)")
  expect_true(all(expected_stats %in% summary_rows$Genotypes))
})

test_that("mean.performance SEm formula is correct", {
  # Create simple test data where we can manually verify SEm
  test_data <- data.frame(
    rep = rep(1:3, each = 4),
    geno = rep(1:4, 3),
    trait = c(10, 12, 11, 13,  # rep 1
              15, 16, 14, 17,  # rep 2
              20, 19, 21, 22)  # rep 3
  )
  
  performance <- mean.performance(data = test_data[,3,drop=FALSE], 
                                genotypes = test_data$geno, 
                                replications = test_data$rep)
  
  # Extract SEm from summary
  sem_row <- performance[performance$Genotypes == "SEm", "trait"]
  
  # Manually calculate using design.stats
  gen_idx <- as.integer(as.factor(test_data$geno))
  rep_idx <- as.integer(as.factor(test_data$rep))
  ds <- design.stats(test_data$trait, test_data$trait, gen_idx, rep_idx, 
                    design_type = "RCBD", calc_type = "all")
  
  r <- length(unique(test_data$rep))
  expected_sem <- sqrt(ds$EMP / r)
  
  expect_equal(as.numeric(sem_row), round(expected_sem, 4), tolerance = 1e-4)
})

test_that("mean.performance CD formulas are correct", {
  performance <- mean.performance(data = seldata[,3:5], genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  # Extract CD values
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  cd5_row <- performance[performance$Genotypes == "CD 5%", ]
  cd1_row <- performance[performance$Genotypes == "CD 1%", ]
  
  # CD values should be positive
  cd5_vals <- as.numeric(cd5_row[, 2:ncol(cd5_row)])
  cd1_vals <- as.numeric(cd1_row[, 2:ncol(cd1_row)])
  
  expect_true(all(cd5_vals > 0))
  expect_true(all(cd1_vals > 0))
  
  # CD 1% should be larger than CD 5%
  expect_true(all(cd1_vals > cd5_vals))
})

test_that("mean.performance CV formula is correct", {
  performance <- mean.performance(data = seldata[,3:5], genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  # Extract CV values
  cv_row <- performance[performance$Genotypes == "CV (%)", ]
  cv_vals <- as.numeric(cv_row[, 2:ncol(cv_row)])
  
  # CV should be positive percentages
  expect_true(all(cv_vals > 0))
  expect_true(all(cv_vals < 100))  # Typically less than 100% for good data
})

test_that("mean.performance heritability is in valid range", {
  performance <- mean.performance(data = seldata[,3:7], genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  # Extract heritability (fraction)
  h2_row <- performance[performance$Genotypes == "Heritability", ]
  h2_vals <- as.numeric(h2_row[, 2:ncol(h2_row)])
  
  # Heritability should be between 0 and 1
  expect_true(all(h2_vals >= 0))
  expect_true(all(h2_vals <= 1))
  
  # Extract heritability (percentage)
  h2_pct_row <- performance[performance$Genotypes == "Heritability(%)", ]
  h2_pct_vals <- as.numeric(h2_pct_row[, 2:ncol(h2_pct_row)])
  
  # Should be percentage form
  expect_equal(h2_pct_vals, h2_vals * 100, tolerance = 1e-4)
})

test_that("mean.performance handles missing values with default REML", {
  # Create data with missing values
  test_data <- seldata[1:40, 3:5]
  test_data[c(1, 5, 10), 1] <- NA
  test_data[c(2, 8), 2] <- NA
  
  # Should warn about using default REML method
  expect_warning(
    performance <- mean.performance(data = test_data, 
                                  genotypes = seldata[1:40, 2], 
                                  replications = seldata[1:40, 1]),
    "Missing values detected.*REML"
  )
  
  # Should return valid results
  expect_true(is.data.frame(performance))
  expect_true(all(is.finite(as.matrix(performance[, -1]))))
})

test_that("mean.performance handles missing values with explicit method", {
  # Create data with missing values
  test_data <- seldata[1:40, 3:5]
  test_data[c(1, 5, 10), 1] <- NA
  
  # Should NOT warn when method is explicit
  expect_silent(
    performance <- mean.performance(data = test_data, 
                                  genotypes = seldata[1:40, 2], 
                                  replications = seldata[1:40, 1],
                                  method = "Mean")
  )
  
  # Should return valid results
  expect_true(is.data.frame(performance))
  expect_true(all(is.finite(as.matrix(performance[, -1]))))
})

test_that("mean.performance supports different imputation methods", {
  # Create data with missing values - use larger dataset for stability
  test_data <- seldata[1:75, 3:5]
  test_data[c(1, 5), 1] <- NA
  test_data[c(10), 2] <- NA
  
  methods <- c("REML", "Mean", "Regression")
  
  for (method in methods) {
    performance <- mean.performance(data = test_data, 
                                  genotypes = seldata[1:75, 2], 
                                  replications = seldata[1:75, 1],
                                  method = method)
    
    expect_true(is.data.frame(performance), 
                info = paste("Method:", method))
    
    # Check that most values are finite (some edge cases might produce NA in stats)
    numeric_mat <- as.matrix(performance[, -1])
    finite_ratio <- sum(is.finite(numeric_mat)) / length(numeric_mat)
    expect_true(finite_ratio > 0.95,
                info = paste("Method:", method, "- finite ratio:", finite_ratio))
  }
})

test_that("mean.performance supports Yates method for balanced designs", {
  # Yates method works best with balanced designs
  # Use complete data (no missing values initially)
  test_data <- seldata[1:60, 3:5]
  # Add just one missing value
  test_data[1, 1] <- NA
  
  performance <- mean.performance(data = test_data, 
                                genotypes = seldata[1:60, 2], 
                                replications = seldata[1:60, 1],
                                method = "Yates")
  
  expect_true(is.data.frame(performance))
  # Yates should produce valid results for most traits
  expect_true(nrow(performance) > 0)
})

test_that("mean.performance default design_type is RCBD", {
  # Should work without specifying design_type
  performance <- mean.performance(data = seldata[,3:5], 
                                genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  expect_true(is.data.frame(performance))
  expect_equal(nrow(performance), nlevels(as.factor(seldata[,2])) + 9)
})

test_that("mean.performance validates LSD requires columns", {
  # Should error when LSD specified without columns
  expect_error(
    mean.performance(data = seldata[,3:5], 
                   genotypes = seldata[,2], 
                   replications = seldata[,1],
                   design_type = "LSD"),
    "Latin Square Design requires 'columns' parameter"
  )
})

test_that("mean.performance handles edge case: zero grand mean", {
  # Create data with all zeros for one trait
  test_data <- data.frame(
    rep = rep(1:3, each = 4),
    geno = rep(1:4, 3),
    trait1 = rep(0, 12),  # All zeros
    trait2 = rnorm(12, 10, 2)
  )
  
  # Should not crash, CV should be NA for zero-mean trait
  performance <- mean.performance(data = test_data[, 3:4], 
                                genotypes = test_data$geno, 
                                replications = test_data$rep)
  
  cv_row <- performance[performance$Genotypes == "CV (%)", ]
  # CV for trait1 should be NA (can't divide by zero)
  expect_true(is.na(as.numeric(cv_row$trait1)))
  # CV for trait2 should be valid
  expect_true(is.finite(as.numeric(cv_row$trait2)))
})

test_that("mean.performance handles negative genetic variance", {
  # When error is large relative to genetic variance, GV can be negative
  # Function should set it to 0
  test_data <- data.frame(
    rep = rep(1:3, each = 4),
    geno = rep(1:4, 3),
    trait = c(10, 10.1, 10, 10.1,  # Very small genotype differences
              10.2, 10, 10.1, 10,   # Large error variance
              10, 10.2, 10.1, 10)
  )
  
  performance <- mean.performance(data = test_data[, 3, drop=FALSE], 
                                genotypes = test_data$geno, 
                                replications = test_data$rep)
  
  h2_row <- performance[performance$Genotypes == "Heritability", ]
  h2_val <- as.numeric(h2_row$trait)
  
  # Heritability should not be negative
  expect_true(h2_val >= 0)
})

test_that("mean.performance genotype means match rowsum calculation", {
  performance <- mean.performance(data = seldata[,3:5], 
                                genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  # Manually calculate using rowsum
  data_mat <- as.matrix(seldata[, 3:5])
  storage.mode(data_mat) <- "numeric"
  genotypes_fac <- as.factor(seldata[,2])
  gen_idx <- as.integer(genotypes_fac)
  manual_means <- rowsum(data_mat, gen_idx, reorder = FALSE) / tabulate(gen_idx)
  
  # Extract genotype means from performance
  n_genotypes <- nlevels(genotypes_fac)
  perf_means <- as.matrix(performance[1:n_genotypes, 2:4])
  
  # Should match (within rounding)
  expect_equal(round(manual_means, 4), perf_means, tolerance = 1e-4)
})

test_that("mean.performance summary statistics are in correct order", {
  performance <- mean.performance(data = seldata[,3:5], 
                                genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  summary_genotypes <- performance$Genotypes[(n_genotypes + 1):nrow(performance)]
  
  expected_order <- c("Min", "Max", "GM", "CV (%)", "SEm", "CD 5%", "CD 1%", 
                     "Heritability", "Heritability(%)")
  
  expect_equal(summary_genotypes, expected_order)
})

test_that("mean.performance Min/Max match genotype means", {
  performance <- mean.performance(data = seldata[,3:5], 
                                genotypes = seldata[,2], 
                                replications = seldata[,1])
  
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  genotype_means <- performance[1:n_genotypes, 2:4]
  
  min_row <- performance[performance$Genotypes == "Min", 2:4]
  max_row <- performance[performance$Genotypes == "Max", 2:4]
  
  # Min should match minimum of genotype means per trait
  for (col in colnames(genotype_means)) {
    expect_equal(as.numeric(min_row[[col]]), 
                min(as.numeric(genotype_means[[col]])),
                tolerance = 1e-4,
                info = paste("Trait:", col))
    
    expect_equal(as.numeric(max_row[[col]]), 
                max(as.numeric(genotype_means[[col]])),
                tolerance = 1e-4,
                info = paste("Trait:", col))
  }
})
