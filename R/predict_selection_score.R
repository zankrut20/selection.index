#' Predict selection index scores
#'
#' @param index_df Data frame returned by lpsi()
#' @param data Raw phenotypic data matrix/data frame (observations x traits)
#' @param genotypes Vector of genotype/treatment labels for each observation
#'
#' @return Data frame of selection index scores by genotype
#' @export
#'
#' @examples
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' cindex <- lpsi(ncomb = 1, pmat = pmat, gmat = gmat, wmat = weight[, -1], wcol = 1)
#' predict_selection_score(cindex, data = seldata[, 3:9], genotypes = seldata[, 2])
#'
predict_selection_score <- function(index_df, data, genotypes) {
  if (!is.data.frame(index_df)) {
    stop("index_df must be a data frame returned by lpsi().")
  }
  if (!"ID" %in% names(index_df)) {
    stop("index_df must contain an ID column.")
  }
  b_cols <- grep("^b\\.", names(index_df), value = TRUE)
  if (length(b_cols) == 0) {
    stop("index_df must contain b.* columns with index coefficients.")
  }

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  if (ncol(data_mat) == 0) {
    stop("data must contain at least one trait column.")
  }

  genotypes_fac <- factor(genotypes, levels = unique(genotypes))
  gen_idx <- as.integer(genotypes_fac)
  if (length(gen_idx) != nrow(data_mat)) {
    stop("genotypes length must match number of rows in data.")
  }

  mean_mat <- cpp_genotype_means(data_mat, gen_idx)
  geno_names <- levels(genotypes_fac)

  n_indices <- nrow(index_df)
  score_mat <- matrix(NA_real_, nrow = length(geno_names), ncol = n_indices)
  score_names <- character(n_indices)

  for (j in seq_len(n_indices)) {
    id_str <- index_df$ID[j]
    idx <- as.integer(trimws(strsplit(id_str, ",")[[1]]))
    if (any(is.na(idx))) {
      stop("ID must contain comma-separated trait indices.")
    }
    if (max(idx) > ncol(mean_mat)) {
      stop("ID indices exceed number of columns in data.")
    }

    b_vals <- as.numeric(index_df[j, b_cols])
    if (length(b_vals) < length(idx)) {
      stop("Number of b coefficients does not match ID length.")
    }
    b_vals <- b_vals[seq_along(idx)] # seq_len(length(idx))]

    trait_means <- mean_mat[, idx, drop = FALSE]
    score_mat[, j] <- rowSums(sweep(trait_means, 2, b_vals, "*"))
    score_names[j] <- paste0("I_", gsub("\\s+", "", gsub(",", "_", id_str)))
  }

  score_names <- make.unique(score_names)

  score_cols <- vector("list", length = n_indices * 2L)
  col_names <- character(n_indices * 2L)
  for (j in seq_len(n_indices)) {
    score_cols[[2L * j - 1L]] <- score_mat[, j]
    col_names[2L * j - 1L] <- score_names[j]
    score_cols[[2L * j]] <- rank(-score_mat[, j], ties.method = "average")
    col_names[2L * j] <- paste0(score_names[j], "_Rank")
  }

  score_df <- data.frame(
    Genotypes = geno_names,
    score_cols,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(score_df) <- c("Genotypes", col_names)

  score_df
}
