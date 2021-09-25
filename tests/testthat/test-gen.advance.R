test_that("Whether gen.advance gives same output", {
  seldata<- rnorm(2,3,4)
  gmat<- matrix(1:16, nrow = 4)
  pmat<- matrix(1:16, nrow = 4)
  weight<- matrix(1.00)
  GA<- gen.advance(phen_mat = pmat[1,1], gen_mat = gmat[1,1], weight_mat = weight)
  result<- structure(2.063, .Dim = c(1L, 1L))
  expect_equal(result, GA)
})
