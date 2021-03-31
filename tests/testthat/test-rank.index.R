test_that("Whether sel.index gives same output", {
  gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
  si<- list()
  si[[1]]<- sel.index(ID = 1, phen_mat = pmat[1,1], gen_mat = gmat[1,1],
                      weight_mat = weight[1,2], GAY = 1.6352)
  si[[2]]<- sel.index(ID = 2, phen_mat = pmat[2,2], gen_mat = gmat[2,2],
                      weight_mat = weight[2,2], GAY = 1.6352)
  si[[3]]<- sel.index(ID = 3, phen_mat = pmat[3,3], gen_mat = gmat[3,3],
                      weight_mat = weight[3,2], GAY = 1.6352)
  si[[4]]<- sel.index(ID = 4, phen_mat = pmat[4,4], gen_mat = gmat[4,4],
                      weight_mat = weight[4,2], GAY = 1.6352)
  rank<- rank.index(list = si, i = 2)

  result<- structure(list(ID = c("2", "1"),
                          b = c(0.5, 0.5),
                          GA = c(1.8221, 1.6352),
                          PRE = c(111.4303, 100.0027),
                          Rank = c(1, 2)),
                     row.names = c("21", "2"),
                     class = "data.frame")

  expect_equal(rank, result)
})
