test_that("gen.varcov and phen.varcov return correct dimensions", {
  gen = gen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  phen = phen_varcov(seldata[,3:9], seldata$treat, seldata$rep)

  expect_equal(nrow(gen), 7)
  expect_equal(ncol(gen), 7)
  expect_equal(nrow(phen), 7)
  expect_equal(ncol(phen), 7)
})

test_that("gen.varcov and phen.varcov return symmetric matrices", {
  gen = gen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  phen = phen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  
  expect_true(isSymmetric(gen))
  expect_true(isSymmetric(phen))
})

test_that("gen.varcov and phen.varcov have appropriate column names", {
  gen = gen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  phen = phen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  
  expect_equal(colnames(gen), colnames(seldata[,3:9]))
  expect_equal(colnames(phen), colnames(seldata[,3:9]))
  expect_equal(rownames(gen), colnames(seldata[,3:9]))
  expect_equal(rownames(phen), colnames(seldata[,3:9]))
})

test_that("gen.varcov returns values less than or equal to phen.varcov", {
  gen = gen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  phen = phen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  
  # Diagonal elements: genotypic variance should be <= phenotypic variance
  expect_true(all(diag(gen) <= diag(phen)))
})

test_that("gen.varcov and phen.varcov handle single trait", {
  gen_single = gen_varcov(seldata[,3, drop=FALSE], seldata$treat, seldata$rep)
  phen_single = phen_varcov(seldata[,3, drop=FALSE], seldata$treat, seldata$rep)
  
  expect_equal(dim(gen_single), c(1, 1))
  expect_equal(dim(phen_single), c(1, 1))
  expect_true(is.finite(gen_single[1,1]))
  expect_true(is.finite(phen_single[1,1]))
})

test_that("gen.varcov and phen.varcov work with different missing value methods", {
  # Create test data with missing values
  test_data <- seldata[,3:9]
  test_data[1, 1] <- NA
  test_data[5, 3] <- NA
  
  # Test all methods
  methods <- c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")
  
  for(method in methods) {
    gen <- gen_varcov(test_data, seldata$treat, seldata$rep, method = method)
    phen <- phen_varcov(test_data, seldata$treat, seldata$rep, method = method)
    
    expect_equal(nrow(gen), 7, info = paste("Method:", method))
    expect_equal(nrow(phen), 7, info = paste("Method:", method))
    expect_true(all(is.finite(gen)), info = paste("Method:", method))
    expect_true(all(is.finite(phen)), info = paste("Method:", method))
  }
})

test_that("gen.varcov and phen.varcov return all finite values", {
  gen = gen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  phen = phen_varcov(seldata[,3:9], seldata$treat, seldata$rep)
  
  expect_true(all(is.finite(gen)))
  expect_true(all(is.finite(phen)))
})
