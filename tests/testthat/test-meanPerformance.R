test_that("Whether result is same as output", {
  performance <- d <- meanPerformance(data = seldata[,3:7], genotypes = seldata[,2], replications = seldata[,1])
  expect_equal(nrow(performance),34)
})
