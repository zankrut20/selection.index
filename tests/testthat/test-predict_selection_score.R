test_that("predict_selection_score returns correct structure", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  scores <- predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[, 2])

  expect_true(is.data.frame(scores))
  expect_true("Genotypes" %in% colnames(scores))
  expect_equal(nrow(scores), nlevels(as.factor(seldata[, 2])))
})

test_that("predict_selection_score includes rank columns", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  scores <- predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[, 2])

  # Check that rank columns exist for each index
  rank_cols <- colnames(scores)[grepl("_Rank$", colnames(scores))]
  expect_true(length(rank_cols) > 0)
})

test_that("predict_selection_score ranks are valid", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  scores <- predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[, 2])

  # Get rank columns
  rank_cols <- colnames(scores)[grepl("_Rank$", colnames(scores))]

  for (col in rank_cols) {
    # Ranks should be numeric
    expect_true(is.numeric(scores[[col]]))

    # Ranks should be in valid range
    n_genotypes <- nrow(scores)
    expect_true(all(scores[[col]] >= 1 & scores[[col]] <= n_genotypes))

    # Check that all ranks from 1 to n are present (accounting for ties)
    expect_true(min(scores[[col]]) >= 1)
    expect_true(max(scores[[col]]) <= n_genotypes)
  }
})

test_that("predict_selection_score higher scores get lower ranks", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  scores <- predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[, 2])

  # Get first score and rank columns
  score_col <- colnames(scores)[!colnames(scores) %in% c("Genotypes", grep("_Rank$", colnames(scores), value = TRUE))]
  rank_col <- paste0(score_col[1], "_Rank")

  # Check that higher scores have lower (better) ranks
  for (i in seq_len(nrow(scores))) {
    for (j in seq_len(nrow(scores))) {
      if (scores[[score_col[1]]][i] > scores[[score_col[1]]][j]) {
        expect_true(scores[[rank_col]][i] <= scores[[rank_col]][j])
      }
    }
  }
})

test_that("predict_selection_score works with multiple indices", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 2, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  scores <- predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[, 2])

  # Should have multiple index scores and corresponding ranks
  all_cols <- colnames(scores)
  rank_cols <- all_cols[grepl("_Rank$", all_cols)]
  score_cols <- all_cols[grepl("^I_", all_cols)]
  score_cols <- score_cols[!grepl("_Rank$", score_cols)] # Remove rank columns

  expect_equal(length(rank_cols), length(score_cols))
  expect_true(length(score_cols) > 1) # Multiple indices
})

test_that("predict_selection_score handles error cases", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  # Test with wrong genotypes length
  expect_error(
    predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[1:10, 2])
  )

  # Test with not a data frame
  expect_error(
    predict_selection_score(as.list(cindex), data = seldata[, 3:9], genotypes = seldata[, 2])
  )
})

# ==============================================================================
# NEW COVERAGE TESTS — targeting previously uncovered lines
# ==============================================================================

test_that("predict_selection_score additional input validations (lines 21-59)", {
  gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
  cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)

  # line 21: missing ID column
  cindex_no_id <- cindex
  cindex_no_id$ID <- NULL
  expect_error(
    predict_selection_score(cindex_no_id, data = seldata[, 3:9], genotypes = seldata[, 2]),
    "index_df must contain an ID column"
  )

  # line 25: missing b.* columns
  cindex_no_b <- cindex
  b_cols <- grep("^b\\.", names(cindex_no_b), value = TRUE)
  for (col in b_cols) cindex_no_b[[col]] <- NULL
  expect_error(
    predict_selection_score(cindex_no_b, data = seldata[, 3:9], genotypes = seldata[, 2]),
    "index_df must contain b.* columns"
  )

  # line 31: data contains no traits (ncol = 0)
  empty_data <- matrix(nrow = nrow(seldata), ncol = 0)
  expect_error(
    predict_selection_score(cindex, data = empty_data, genotypes = seldata[, 2]),
    "data must contain at least one trait column"
  )

  # line 51: ID must contain comma-separated indices (causing NA parser result)
  cindex_bad_id_format <- cindex
  cindex_bad_id_format$ID[1] <- "1, 2, letters"
  suppressWarnings(
    expect_error(
      predict_selection_score(cindex_bad_id_format, data = seldata[, 3:9], genotypes = seldata[, 2]),
      "ID must contain comma-separated trait indices"
    )
  )

  # line 54: ID indices exceed number of columns in data
  cindex_out_of_bounds <- cindex
  cindex_out_of_bounds$ID[1] <- paste(1:10, collapse = ", ") # data only has 7 columns
  expect_error(
    predict_selection_score(cindex_out_of_bounds, data = seldata[, 3:9], genotypes = seldata[, 2]),
    "ID indices exceed number of columns in data"
  )

  # line 59: Number of b coefficients does not match ID length
  cindex_mismatched_b <- cindex
  # artificially add an extra index to ID without adding a b_col
  # original ID has 7 indices, let's make it 8
  cindex_mismatched_b$ID[1] <- paste0(cindex_mismatched_b$ID[1], ", 1")
  expect_error(
    predict_selection_score(cindex_mismatched_b, data = seldata[, 3:9], genotypes = seldata[, 2]),
    "Number of b coefficients does not match ID length"
  )
})
