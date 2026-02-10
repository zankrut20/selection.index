test_that("rlpsi returns expected structure and satisfies constraints", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]

  # Test with auto-created C (user-friendly way)
  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = 1)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))
  expect_true("C" %in% names(result))

  summary_df <- result$summary
  expect_true(is.data.frame(summary_df))
  expect_true(all(c("GA", "PRE", "Delta_G", "rHI", "hI2") %in% colnames(summary_df)))

  # Check b is a clean numeric vector
  expect_true(is.numeric(result$b))
  expect_false(is.matrix(result$b))
  expect_equal(length(result$b), ncol(pmat))
  expect_equal(length(result$Delta_G), ncol(pmat))

  # Constraint satisfaction: C' Delta_G should be close to zero
  constraint_value <- as.numeric(t(result$C) %*% as.numeric(result$Delta_G))
  expect_true(abs(constraint_value) < 1e-4)
})

test_that("rlpsi works with custom C matrix (backward compatibility)", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  C <- diag(ncol(pmat))[, 1, drop = FALSE]

  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, C = C)

  expect_true(is.list(result))
  expect_true("C" %in% names(result))
  
  # Verify constraint is satisfied
  constraint_value <- as.numeric(t(result$C) %*% as.numeric(result$Delta_G))
  expect_true(abs(constraint_value) < 1e-4)
})

test_that("rlpsi can restrict multiple traits", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]

  # Restrict traits 1 and 3
  result <- selection.index:::rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = c(1, 3))

  expect_equal(ncol(result$C), 2)
  
  # Both constraints should be satisfied
  constraint_values <- as.numeric(t(result$C) %*% as.numeric(result$Delta_G))
  expect_true(all(abs(constraint_values) < 1e-4))
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

  # Check b is a clean numeric vector
  expect_true(is.numeric(result$b))
  expect_false(is.matrix(result$b))
  expect_equal(length(result$b), ncol(pmat))

  ratios <- as.numeric(result$Delta_G) / k
  ratios <- ratios[is.finite(ratios)]
  expect_true((max(ratios) - min(ratios)) < 1e-4)
})

test_that("dg_lpsi returns expected structure and achieves desired gains", {
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  d <- seq_len(ncol(pmat))

  # Note: DG-LPSI doesn't use economic weights (wmat)
  result <- selection.index:::dg_lpsi(pmat, gmat, d)

  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("b" %in% names(result))
  expect_true("Delta_G" %in% names(result))
  expect_true("desired_gains" %in% names(result))

  summary_df <- result$summary
  
  # DG-LPSI should NOT have GA and PRE (no economic weights)
  expect_false("GA" %in% colnames(summary_df))
  expect_false("PRE" %in% colnames(summary_df))
  
  # DG-LPSI should have hI2 and rHI
  expect_true("hI2" %in% colnames(summary_df))
  expect_true("rHI" %in% colnames(summary_df))
  expect_true("Delta_G" %in% colnames(summary_df))

  # Check b is a clean numeric vector
  expect_true(is.numeric(result$b))
  expect_false(is.matrix(result$b))
  expect_equal(length(result$b), ncol(pmat))

  # Check that desired gains match input
  expect_equal(as.numeric(result$desired_gains), d)
  
  # Check that realized gains are proportional to desired gains
  ratios <- as.numeric(result$Delta_G) / d
  ratios <- ratios[is.finite(ratios)]
  expect_true((max(ratios) - min(ratios)) < 1e-4)
})

test_that("ppg_lpsi detects singular matrices", {
  # Create a rank-deficient G matrix (last trait is linear combination of first two)
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  wmat <- weight[,-1]
  
  # Make G rank-deficient by setting last column/row as linear combination
  gmat_singular <- gmat
  gmat_singular[, ncol(gmat)] <- gmat_singular[, 1] + 0.5 * gmat_singular[, 2]
  gmat_singular[ncol(gmat), ] <- gmat_singular[1, ] + 0.5 * gmat_singular[2, ]
  
  k <- rep(1, ncol(pmat))
  
  # Should detect singular matrix and throw error
  expect_error(
    selection.index:::ppg_lpsi(pmat, gmat_singular, k, wmat = wmat, wcol = 1),
    "Singular matrix detected"
  )
})

test_that("dg_lpsi detects singular matrices", {
  # Create a rank-deficient G matrix
  gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
  
  # Make G rank-deficient
  gmat_singular <- gmat
  gmat_singular[, ncol(gmat)] <- gmat_singular[, 1] + 0.5 * gmat_singular[, 2]
  gmat_singular[ncol(gmat), ] <- gmat_singular[1, ] + 0.5 * gmat_singular[2, ]
  
  d <- seq_len(ncol(pmat))
  
  # Should detect singular matrix and throw error
  expect_error(
    selection.index:::dg_lpsi(pmat, gmat_singular, d),
    "Singular matrix detected"
  )
})
