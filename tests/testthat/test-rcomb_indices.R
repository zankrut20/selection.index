test_that("rcomb.indices returns correct dimensions", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex), 6)
})

test_that("rcomb.indices returns data frame with correct columns", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  expect_true(is.data.frame(cindex))
  expect_true("ID" %in% colnames(cindex))
  expect_true(any(grepl("^b\\.", colnames(cindex))))
  expect_true("GA" %in% colnames(cindex))
})

test_that("rcomb.indices correctly removes specified trait", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Remove trait 1, ncomb = 1 should give 6 combinations (7 - 1)
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex), 6)
  
  # Check that none of the IDs contain "1"
  expect_true(all(!grepl("^1$", cindex$ID)))
})

test_that("rcomb.indices works with different ncomb values", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # ncomb = 2, remove trait 1: choose(6, 2) = 15 combinations
  cindex2<- rcomb_indices(ncomb = 2, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex2), 15)
  
  # ncomb = 3, remove trait 1: choose(6, 3) = 20 combinations
  cindex3<- rcomb_indices(ncomb = 3, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex3), 20)
})

test_that("rcomb.indices calculates PRE when GAY is provided", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1, GAY = 1.075)
  
  expect_true("PRE" %in% colnames(cindex))
  expect_true(all(is.numeric(cindex$PRE)))
  expect_true(all(is.finite(cindex$PRE)))
})

test_that("rcomb.indices returns numeric GA values", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  expect_true(all(is.numeric(cindex$GA)))
  expect_true(all(is.finite(cindex$GA)))
})

test_that("rcomb.indices works with different removed traits", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Remove trait 1
  cindex1<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex1), 6)
  
  # Remove trait 3
  cindex3<- rcomb_indices(ncomb = 1, i = 3, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex3), 6)
  
  # Check that trait 3 is not in IDs
  expect_true(all(!grepl("^3$", cindex3$ID)))
})

test_that("rcomb.indices works with different weight columns", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Test with wcol = 1
  cindex1<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex1), 6)
  
  # Test with wcol = 2 (if available)
  if(ncol(weight[,-1]) >= 2) {
    cindex2<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 2)
    expect_equal(nrow(cindex2), 6)
  }
})
test_that("rcomb.indices calculates new metrics: Delta_G, rHI, hI2", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  # Check that new columns exist
  expect_true("Delta_G" %in% colnames(cindex))
  expect_true("rHI" %in% colnames(cindex))
  expect_true("hI2" %in% colnames(cindex))
})

test_that("rcomb.indices new metrics are numeric and valid", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  # Check Delta_G values
  expect_true(all(is.numeric(cindex$Delta_G)))
  expect_true(all(is.finite(cindex$Delta_G)))
  expect_true(all(cindex$Delta_G >= 0))  # Delta_G should be non-negative
  
  # Check rHI values (should be between 0 and 1)
  expect_true(all(is.numeric(cindex$rHI)))
  expect_true(all(is.finite(cindex$rHI)))
  expect_true(all(cindex$rHI >= 0 & cindex$rHI <= 1))
  
  # Check hI2 values (should be between 0 and 1)
  expect_true(all(is.numeric(cindex$hI2)))
  expect_true(all(is.finite(cindex$hI2)))
  expect_true(all(cindex$hI2 >= 0 & cindex$hI2 <= 1))
})

test_that("rcomb.indices rHI equals sqrt(hI2)", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb_indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  # rHI should equal sqrt(hI2) mathematically
  expected_rHI <- sqrt(cindex$hI2)
  expect_equal(cindex$rHI, expected_rHI, tolerance = 1e-4)
})