test_that("Whether both functions gives us the same output", {
  gen = gen.varcov(seldata[,3:9], seldata$treat, seldata$rep)
  phen = phen.varcov(seldata[,3:9], seldata$treat, seldata$rep)

  expect_equal(nrow(gen), 7)
  expect_equal(nrow(phen), 7)
})
