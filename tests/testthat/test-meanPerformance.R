test_that("meanPerformance returns correct dimensions", {
  performance <- meanPerformance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  expect_equal(nrow(performance), 34)
})

test_that("meanPerformance returns data frame with correct columns", {
  performance <- meanPerformance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  expect_true(is.data.frame(performance))
  expect_true("Genotypes" %in% colnames(performance))
  # Column names should match trait names from input data
  expect_true(all(colnames(seldata[,3:7]) %in% colnames(performance)))
})

test_that("meanPerformance handles all traits", {
  performance <- meanPerformance(data = seldata[,3:9], genotypes = seldata[,2], replications = seldata[,1])
  
  # Should have rows for genotypes + 9 summary statistics (Min, Max, GM, CV, SEm, CD 5%, CD 1%, Heritability, Heritability%)
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  expect_equal(nrow(performance), n_genotypes + 9)
  
  # Check that trait columns exist
  expect_true(all(colnames(seldata[,3:9]) %in% colnames(performance)))
})

test_that("meanPerformance genotype means are correct", {
  performance <- meanPerformance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  # Manually calculate mean for first genotype and first trait
  first_geno <- levels(as.factor(seldata[,2]))[1]
  first_trait <- colnames(seldata)[3]
  manual_mean <- mean(seldata[seldata[,2] == first_geno, 3], na.rm = TRUE)
  
  perf_mean <- performance[performance$Genotypes == first_geno, first_trait]
  
  # Should be close (allowing for rounding)
  expect_true(abs(manual_mean - as.numeric(perf_mean)) < 0.01)
})

test_that("meanPerformance works with single trait", {
  performance <- meanPerformance(data = seldata[,3, drop=FALSE], genotypes = seldata[,2], replications = seldata[,1])
  
  expect_true(is.data.frame(performance))
  # Should have genotypes + 9 summary statistics rows
  n_genotypes <- nlevels(as.factor(seldata[,2]))
  expect_equal(nrow(performance), n_genotypes + 9)
})

test_that("meanPerformance includes all genotypes", {
  performance <- meanPerformance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  unique_genos <- unique(seldata[,2])
  # First n rows should be genotype means, last 9 rows are summary statistics
  n_genotypes <- length(unique_genos)
  genotype_rows <- performance$Genotypes[1:n_genotypes]
  
  # Check that all unique genotypes appear in the first n rows
  expect_true(all(unique_genos %in% genotype_rows))
})

test_that("meanPerformance calculates statistics correctly", {
  performance <- meanPerformance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  
  # Check that all numeric columns contain finite values
  numeric_cols <- sapply(performance, is.numeric)
  numeric_data <- performance[, numeric_cols]
  expect_true(all(is.finite(as.matrix(numeric_data))))
})
