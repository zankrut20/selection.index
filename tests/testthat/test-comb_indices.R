test_that("comb.indices returns correct dimensions", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex), 7)
})

test_that("comb.indices returns data frame with correct columns", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  expect_true(is.data.frame(cindex))
  expect_true("ID" %in% colnames(cindex))
  expect_true(any(grepl("^b\\.", colnames(cindex))))
  expect_true("GA" %in% colnames(cindex))
})

test_that("comb.indices handles different ncomb values correctly", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Test single character combination
  cindex1<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex1), 7)  # choose(7, 1) = 7
  
  # Test pair combinations
  cindex2<- comb_indices(ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex2), 21)  # choose(7, 2) = 21
})

test_that("comb.indices calculates PRE when GAY is provided", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1, GAY = 1.075)
  
  expect_true("PRE" %in% colnames(cindex))
  expect_true(all(is.numeric(cindex$PRE)))
})

test_that("comb.indices works with different weight columns", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Test with wcol = 1
  cindex1<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex1), 7)
  
  # Test with wcol = 2 (if available)
  if(ncol(weight[,-1]) >= 2) {
    cindex2<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 2)
    expect_equal(nrow(cindex2), 7)
  }
})

test_that("comb.indices returns numeric GA values", {
  gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  
  expect_true(all(is.numeric(cindex$GA)))
  expect_true(all(is.finite(cindex$GA)))
})
