test_that("Whether result is same as output", {
  sr<- sel.score.rank(data = seldata[1:9,3],
                      bmat = 0.6455, genotype = seldata[1:9,2])
  expect_equal(nrow(sr),3)
})
