test_that("rlpsi returns expected structure and satisfies constraints", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  C <- diag(ncol(pmat))[, 1, drop = FALSE]

  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, C = C)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))

  summary_df <- result$summary
  expect_true(is.data.frame(summary_df))
  expect_true(all(c("GA", "PRE", "Delta_G", "rHI", "hI2") %in% colnames(summary_df)))

  expect_equal(length(result$b), ncol(pmat))
  expect_equal(length(result$Delta_G), ncol(pmat))

  # Constraint satisfaction: C' Delta_G should be close to zero
  constraint_value <- as.numeric(t(C) %*% as.numeric(result$Delta_G))
  expect_true(abs(constraint_value) < 1e-4)
})

test_that("ppg_lpsi returns expected structure and proportional gains", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  k <- rep(1, ncol(pmat))

  result <- selection.index:::ppg_lpsi(pmat, gmat, k, wmat = wmat, wcol = 1)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))
  expect_true("phi" %in% names(result))

  summary_df <- result$summary
  expect_true("phi" %in% colnames(summary_df))
  expect_true(is.numeric(summary_df$phi))

  ratios <- as.numeric(result$Delta_G) / k
  ratios <- ratios[is.finite(ratios)]
  expect_true((max(ratios) - min(ratios)) < 1e-4)
})

test_that("dg_lpsi returns expected structure and proportional desired gains", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  d <- seq_len(ncol(pmat))

  result <- selection.index:::dg_lpsi(pmat, gmat, d, wmat = wmat, wcol = 1)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))

  ratios <- as.numeric(result$Delta_G) / d
  ratios <- ratios[is.finite(ratios)]
  expect_true((max(ratios) - min(ratios)) < 1e-4)
})
