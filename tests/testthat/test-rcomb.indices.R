test_that("Whether result is same as output", {
  gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- rcomb.indices(ncomb = 1, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,-1], wcol = 1)
  expect_equal(nrow(cindex),6)
})
