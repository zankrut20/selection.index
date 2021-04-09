test_that("Whether result is same as output", {
  gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
  cindex<- comb.indices(ncomb = 1, phen_mat = pmat, gen_mat = gmat, weight_mat = weight[,-1], weight_col = 1)
  expect_equal(nrow(cindex),7)
})
