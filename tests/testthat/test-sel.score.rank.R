test_that("Whether result is same as output", {
  sr<- sel.score.rank(data = seldata[1:9,3],
                      bmat = 0.6455, genotype = seldata[1:9,2])
  result<- structure(
    list(
      Genotype = c("G1", "G2", "G3"),
      Selection.score = c(3.53232661666667, 4.31746978333333,
                          2.93848813333333),
      Rank = c(2, 1, 3)), class = "data.frame", row.names = c(NA, -3L))
  expect_equal(result, sr)
})
