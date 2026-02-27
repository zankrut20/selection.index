#' Missing Value Estimation for RCBD, Latin Square, and Split Plot Designs
#'
#' Estimates missing values in Randomized Complete Block Design (RCBD),
#' Latin Square Design (LSD), or Split Plot Design (SPD) using one of six methods:
#' REML, Yates, Healy, Regression, Mean, or Bartlett.
#'
#' @param data_mat Numeric matrix with observations (rows) by traits (columns).
#'   May contain missing values (NA, NaN, Inf).
#' @param gen_idx Integer vector indicating genotype/treatment index for each observation (sub-plot in SPD).
#' @param rep_idx Integer vector indicating replicate/block (RCBD) or row (LSD) index for each observation.
#' @param col_idx Integer vector indicating column index for each observation (required for LSD only).
#' @param main_idx Integer vector indicating main plot treatment index for each observation (required for SPD only).
#' @param design_type Character string specifying experimental design: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot).
#' @param method Character string specifying the estimation method. One of:
#'   \describe{
#'     \item{REML}{Restricted Maximum Likelihood with variance components and BLUP
#'       (Best Linear Unbiased Prediction). Most robust for complex missing patterns.
#'       Default method. Iterations: 100. (RCBD, LSD only)}
#'     \item{Yates}{Traditional iterative formula: (r*T + t*B - G) / ((r-1)(t-1)).
#'       Simple and fast. Good for simple missing patterns. Iterations: 50. (RCBD, LSD only)}
#'     \item{Healy}{Healy & Westmacott method using weighted adjustment of treatment
#'       and block means. More stable than Yates for multiple missing values. Iterations: 50. (RCBD, LSD only)}
#'     \item{Regression}{Linear regression on treatment and block effects with QR
#'       decomposition. Non-iterative, single-pass estimation. Fast and stable. (RCBD, LSD only)}
#'     \item{Mean}{Simple mean substitution using treatment and block effects.
#'       Non-iterative, fastest method. Good for quick estimation. (RCBD, LSD, SPD)}
#'     \item{Bartlett}{ANCOVA using other traits as covariates with QR decomposition.
#'       Best when traits are correlated. Iterations: 30. (RCBD, LSD only)}
#'   }
#' @param tolerance Numeric convergence criterion. Iteration stops when maximum
#'   change in estimated values is below this threshold. Default: 1e-6.
#'
#' @return A numeric matrix of the same dimensions as \code{data_mat} with all
#'   missing values replaced by estimates.
#'
#' @details
#' The function handles missing values in RCBD, LSD, or SPD experiments by iteratively
#' estimating them until convergence. For RCBD, uses 2-way blocking (genotypes × blocks);
#' for LSD, uses 3-way blocking (genotypes × rows × columns); for SPD, uses nested
#' structure (blocks > main plots > sub-plots). Each method has different strengths:
#'
#' \strong{REML Method:}
#' Uses variance component estimation with restricted maximum likelihood.
#' Applies shrinkage to predictions based on estimated variance ratios.
#' Most computationally intensive but handles complex patterns well.
#'
#' \strong{Yates Method:}
#' Uses the classical formula accounting for treatment, block, and grand totals.
#' Fast and simple, suitable for balanced designs with few missing values.
#'
#' \strong{Healy & Westmacott Method:}
#' Uses weighted adjustment based on residuals from treatment and block effects.
#' More stable than Yates when multiple values are missing from same treatment or block.
#' Particularly effective for unbalanced missing patterns.
#'
#' \strong{Regression Method:}
#' Fits linear model using complete observations with treatment and block as factors.
#' Non-iterative single-pass estimation using QR decomposition for numerical stability.
#' Fast and deterministic - no convergence iterations needed.
#'
#' \strong{Mean Substitution Method:}
#' Replaces missing values with treatment mean + block mean - grand mean.
#' Non-iterative, fastest method. Simple additive model without regression.
#' Best for quick estimation when precision is less critical.
#'
#' \strong{Bartlett Method:}
#' Performs ANCOVA using other traits as covariates to predict missing values.
#' Leverages trait correlations for better predictions when traits are related.
#' Uses QR decomposition for numerical stability.
#'
#' \strong{Method Availability by Design:}
#' \itemize{
#'   \item RCBD: REML, Yates, Healy, Regression, Mean, Bartlett
#'   \item LSD:  REML, Yates, Healy, Regression, Mean, Bartlett
#'   \item SPD:  Mean only (other methods fall back to Mean)
#' }
#'
#' @references
#' Yates, F. (1933). The analysis of replicated experiments when the field results
#' are incomplete. \emph{Empire Journal of Experimental Agriculture}, 1, 129-142.
#'
#' Healy, M. J. R., & Westmacott, M. (1956). Missing values in experiments analysed
#' on automatic computers. \emph{Applied Statistics}, 5(3), 203-206.
#'
#' Bartlett, M. S. (1937). Some examples of statistical methods of research in
#' agriculture and applied biology. \emph{Supplement to the Journal of the Royal
#' Statistical Society}, 4(2), 137-183.
#'
#' @keywords internal
#' @noRd
missing_value_estimation <- function(data_mat, gen_idx, rep_idx, col_idx = NULL, main_idx = NULL,
                                     design_type = c("RCBD", "LSD", "SPD"),
                                     method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"),
                                     tolerance = 1e-6) {
  design_type <- match.arg(design_type)
  method <- match.arg(method)

  if (design_type == "LSD" && is.null(col_idx)) {
    stop("Latin Square Design requires 'col_idx' parameter")
  }

  if (design_type == "SPD" && is.null(main_idx)) {
    stop("Split Plot Design requires 'main_idx' parameter")
  }

  if (design_type == "SPD" && method != "Mean") {
    warning("For Split Plot Design, switching to 'Mean' method for missing value estimation. ",
      "SPD has nested structure best handled by mean substitution.",
      call. = FALSE
    )
    method <- "Mean"
  }


  ncols <- ncol(data_mat)
  genotype <- length(unique(gen_idx))
  repli <- length(unique(rep_idx))
  ncol_blocks <- if (!is.null(col_idx)) length(unique(col_idx)) else 0


  if (!any(!is.finite(data_mat))) {
    return(data_mat) # No missing values, return as is
  }

  data_imputed <- data_mat

  max_iter <- if (method == "REML") {
    100
  } else if (method == "Yates" || method == "Healy") {
    50
  } else if (method == "Regression" || method == "Mean") {
    1 # Non-iterative
  } else {
    30
  }

  if (method == "Yates") {

    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))

      if (length(missing_idx) > 0) {
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean

        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]

          gen_sums <- grouped_sums(matrix(data_imputed[, col], ncol = 1), gen_idx)[, 1]
          rep_sums <- grouped_sums(matrix(data_imputed[, col], ncol = 1), rep_idx)[, 1]
          G_total <- sum(data_imputed[, col], na.rm = TRUE)

          if (design_type == "LSD") {
            col_sums <- grouped_sums(matrix(data_imputed[, col], ncol = 1), col_idx)[, 1]
          }

          for (idx in missing_idx) {
            g <- gen_idx[idx]
            r <- rep_idx[idx]

            current_val <- data_imputed[idx, col]
            T_total <- gen_sums[g] - current_val
            R_total <- rep_sums[r] - current_val

            if (design_type == "RCBD") {
              numerator <- repli * T_total + genotype * R_total - (G_total - current_val)
              denominator <- (repli - 1) * (genotype - 1)
            } else {
              c <- col_idx[idx]
              C_total <- col_sums[c] - current_val
              numerator <- genotype * T_total + genotype * R_total + genotype * C_total - 2 * (G_total - current_val)
              denominator <- (genotype - 1) * (genotype - 2)
            }

            data_imputed[idx, col] <- numerator / denominator
          }

          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
  } else if (method == "Healy") {

    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))

      if (length(missing_idx) > 0) {
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean

        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]

          gen_sums <- grouped_sums(matrix(data_imputed[, col], ncol = 1), gen_idx)[, 1]
          rep_sums <- grouped_sums(matrix(data_imputed[, col], ncol = 1), rep_idx)[, 1]

          complete_mask <- matrix(as.numeric(is.finite(data_mat[, col])), ncol = 1)
          gen_counts <- grouped_sums(complete_mask, gen_idx)[, 1]
          rep_counts <- grouped_sums(complete_mask, rep_idx)[, 1]

          treat_means <- gen_sums / pmax(gen_counts, 1)
          block_means <- rep_sums / pmax(rep_counts, 1)
          grand_mean <- mean(data_imputed[, col])

          if (design_type == "LSD") {
            col_sums <- grouped_sums(matrix(data_imputed[, col], ncol = 1), col_idx)[, 1]
            col_counts <- grouped_sums(complete_mask, col_idx)[, 1]
            col_means <- col_sums / pmax(col_counts, 1)
          }

          for (idx in missing_idx) {
            g <- gen_idx[idx]
            r <- rep_idx[idx]

            n_miss_treat <- sum(!is.finite(data_mat[gen_idx == g, col]))
            n_miss_block <- sum(!is.finite(data_mat[rep_idx == r, col]))

            if (design_type == "RCBD") {
              w_treat <- (repli - n_miss_treat) / repli
              w_block <- (genotype - n_miss_block) / genotype

              w_sum <- w_treat + w_block
              if (w_sum > 0) {
                w_treat <- w_treat / w_sum
                w_block <- w_block / w_sum
              } else {
                w_treat <- 0.5
                w_block <- 0.5
              }

              treat_effect <- treat_means[g] - grand_mean
              block_effect <- block_means[r] - grand_mean

              data_imputed[idx, col] <- grand_mean +
                w_treat * treat_effect +
                w_block * block_effect
            } else {
              c <- col_idx[idx]
              n_miss_col <- sum(!is.finite(data_mat[col_idx == c, col]))

              w_treat <- (genotype - n_miss_treat) / genotype
              w_row <- (genotype - n_miss_block) / genotype
              w_col <- (genotype - n_miss_col) / genotype

              w_sum <- w_treat + w_row + w_col
              if (w_sum > 0) {
                w_treat <- w_treat / w_sum
                w_row <- w_row / w_sum
                w_col <- w_col / w_sum
              } else {
                w_treat <- 1 / 3
                w_row <- 1 / 3
                w_col <- 1 / 3
              }

              treat_effect <- treat_means[g] - grand_mean
              row_effect <- block_means[r] - grand_mean
              col_effect <- col_means[c] - grand_mean

              data_imputed[idx, col] <- grand_mean +
                w_treat * treat_effect +
                w_row * row_effect +
                w_col * col_effect
            }
          }

          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
  } else if (method == "Regression") {



    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))

      if (length(missing_idx) > 0) {
        complete_idx <- which(is.finite(data_mat[, col]))

        min_obs <- if (design_type == "RCBD") {
          genotype + repli
        } else {
          genotype + repli + ncol_blocks
        }

        if (length(complete_idx) > min_obs) {
          n_complete <- length(complete_idx)
          n_missing <- length(missing_idx)

          y_complete <- data_mat[complete_idx, col]
          treat_complete <- gen_idx[complete_idx]
          block_complete <- rep_idx[complete_idx]

          X_treat <- matrix(0, n_complete, genotype - 1)
          for (g in seq_len(genotype - 1)) {
            X_treat[treat_complete == g, g] <- 1
          }
          X_treat[treat_complete == genotype, ] <- -1

          X_block <- matrix(0, n_complete, repli - 1)
          for (r in seq_len(repli - 1)) {
            X_block[block_complete == r, r] <- 1
          }
          X_block[block_complete == repli, ] <- -1

          if (design_type == "LSD") {
            col_complete <- col_idx[complete_idx]
            X_col <- matrix(0, n_complete, ncol_blocks - 1)
            for (c in seq_len(ncol_blocks - 1)) {
              X_col[col_complete == c, c] <- 1
            }
            X_col[col_complete == ncol_blocks, ] <- -1
            X <- cbind(1, X_treat, X_block, X_col)
          } else {
            X <- cbind(1, X_treat, X_block)
          }

          qr_fit <- qr(X)

          if (qr_fit$rank == ncol(X)) {
            beta <- qr.coef(qr_fit, y_complete)

            X_missing <- matrix(0, n_missing, ncol(X))
            X_missing[, 1] <- 1 # Intercept column

            treat_missing <- gen_idx[missing_idx]
            for (g in seq_len(genotype - 1)) {
              X_missing[treat_missing == g, 1 + g] <- 1
            }
            X_missing[treat_missing == genotype, 2:genotype] <- -1

            block_missing <- rep_idx[missing_idx]
            block_start <- genotype
            for (r in seq_len(repli - 1)) {
              X_missing[block_missing == r, block_start + r] <- 1
            }
            X_missing[block_missing == repli, (block_start + 1):(block_start + repli - 1)] <- -1

            if (design_type == "LSD") {
              col_missing <- col_idx[missing_idx]
              col_start <- genotype + repli - 1
              for (c in seq_len(ncol_blocks - 1)) {
                X_missing[col_missing == c, col_start + c] <- 1
              }
              X_missing[col_missing == ncol_blocks, (col_start + 1):(col_start + ncol_blocks - 1)] <- -1
            }

            data_imputed[missing_idx, col] <- X_missing %*% beta
          } else {
            complete_data <- matrix(data_mat[complete_idx, col], ncol = 1)
            treat_sums <- grouped_sums(complete_data, gen_idx[complete_idx])[, 1]
            treat_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), gen_idx[complete_idx])[, 1]
            treat_means <- treat_sums / pmax(treat_counts, 1)

            block_sums <- grouped_sums(complete_data, rep_idx[complete_idx])[, 1]
            block_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), rep_idx[complete_idx])[, 1]
            block_means <- block_sums / pmax(block_counts, 1)

            grand_mean <- cpp_grand_means(complete_data)[1]

            if (design_type == "LSD") {
              col_sums <- grouped_sums(complete_data, col_idx[complete_idx])[, 1]
              col_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), col_idx[complete_idx])[, 1]
              col_means <- col_sums / pmax(col_counts, 1)

              for (idx in missing_idx) {
                g <- gen_idx[idx]
                r <- rep_idx[idx]
                c <- col_idx[idx]
                treat_eff <- if (g <= length(treat_means)) treat_means[g] - grand_mean else 0
                row_eff <- if (r <= length(block_means)) block_means[r] - grand_mean else 0
                col_eff <- if (c <= length(col_means)) col_means[c] - grand_mean else 0
                data_imputed[idx, col] <- grand_mean + treat_eff + row_eff + col_eff
              }
            } else {
              for (idx in missing_idx) {
                g <- gen_idx[idx]
                r <- rep_idx[idx]
                treat_eff <- if (g <= length(treat_means)) treat_means[g] - grand_mean else 0
                block_eff <- if (r <= length(block_means)) block_means[r] - grand_mean else 0
                data_imputed[idx, col] <- grand_mean + treat_eff + block_eff
              }
            }
          }
        } else {
          data_imputed[missing_idx, col] <- mean(data_mat[, col], na.rm = TRUE)
        }
      }
    }
  } else if (method == "Mean") {



    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))

      if (length(missing_idx) > 0) {
        complete_idx <- which(is.finite(data_mat[, col]))

        if (length(complete_idx) > 0) {
          y_complete <- data_mat[complete_idx, col]
          treat_complete <- gen_idx[complete_idx]
          block_complete <- rep_idx[complete_idx]

          complete_data <- matrix(y_complete, ncol = 1)

          treat_sums <- grouped_sums(complete_data, treat_complete)[, 1]
          gen_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), treat_complete)[, 1]
          treat_means <- treat_sums / pmax(gen_counts, 1)

          block_sums <- grouped_sums(complete_data, block_complete)[, 1]
          rep_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), block_complete)[, 1]
          block_means <- block_sums / pmax(rep_counts, 1)

          grand_mean <- cpp_grand_means(complete_data)[1]

          if (design_type == "LSD") {
            col_complete <- col_idx[complete_idx]
            col_sums <- grouped_sums(complete_data, col_complete)[, 1]
            col_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), col_complete)[, 1]
            col_means <- col_sums / pmax(col_counts, 1)

            for (idx in missing_idx) {
              g <- gen_idx[idx]
              r <- rep_idx[idx]
              c <- col_idx[idx]

              treat_effect <- if (g <= length(treat_means)) treat_means[g] - grand_mean else 0

              row_effect <- if (r <= length(block_means)) block_means[r] - grand_mean else 0

              col_effect <- if (c <= length(col_means)) col_means[c] - grand_mean else 0

              data_imputed[idx, col] <- grand_mean + treat_effect + row_effect + col_effect
            }
          } else if (design_type == "SPD") {
            main_complete <- main_idx[complete_idx]
            main_sums <- grouped_sums(complete_data, main_complete)[, 1]
            main_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), main_complete)[, 1]
            main_means <- main_sums / pmax(main_counts, 1)

            for (idx in missing_idx) {
              g <- gen_idx[idx] # Sub-plot treatment (genotype)
              r <- rep_idx[idx] # Block/replication
              m <- main_idx[idx] # Main plot treatment

              block_effect <- if (r <= length(block_means)) block_means[r] - grand_mean else 0

              main_effect <- if (m <= length(main_means)) main_means[m] - grand_mean else 0

              treat_effect <- if (g <= length(treat_means)) treat_means[g] - grand_mean else 0

              data_imputed[idx, col] <- grand_mean + block_effect + main_effect + treat_effect
            }
          } else {
            for (idx in missing_idx) {
              g <- gen_idx[idx]
              r <- rep_idx[idx]

              treat_effect <- if (g <= length(treat_means)) treat_means[g] - grand_mean else 0

              block_effect <- if (r <= length(block_means)) block_means[r] - grand_mean else 0

              data_imputed[idx, col] <- grand_mean + treat_effect + block_effect
            }
          }
        } else {
          data_imputed[missing_idx, col] <- 0
        }
      }
    }
  } else if (method == "REML") {

    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))

      if (length(missing_idx) > 0) {
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean

        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]

          complete_idx <- which(is.finite(data_mat[, col]))
          if (length(complete_idx) > 3) {
            y_complete <- data_imputed[complete_idx, col]
            treat_complete <- gen_idx[complete_idx]
            block_complete <- rep_idx[complete_idx]

            complete_data <- matrix(y_complete, ncol = 1)
            grand_mean <- cpp_grand_means(complete_data)[1]

            treat_sums <- grouped_sums(complete_data, treat_complete)[, 1]
            treat_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), treat_complete)[, 1]
            unique_treats <- sort(unique(treat_complete))
            treat_means_vec <- treat_sums / pmax(treat_counts, 1)
            treat_means <- numeric(genotype)
            treat_means[unique_treats] <- treat_means_vec

            block_sums <- grouped_sums(complete_data, block_complete)[, 1]
            block_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), block_complete)[, 1]
            unique_blocks <- sort(unique(block_complete))
            block_means_vec <- block_sums / pmax(block_counts, 1)
            block_means <- numeric(repli)
            block_means[unique_blocks] <- block_means_vec

            ss_treat <- sum(treat_counts * (treat_means_vec - grand_mean)^2)
            ss_block <- sum(block_counts * (block_means_vec - grand_mean)^2)

            if (design_type == "RCBD") {
              residuals <- y_complete - treat_means[treat_complete] -
                block_means[block_complete] + grand_mean
              ss_error <- sum(residuals^2)

              df_treat <- genotype - 1
              df_block <- repli - 1
              df_error <- length(complete_idx) - genotype - repli + 1

              if (df_error > 0) {
                ms_treat <- ss_treat / df_treat
                ms_block <- ss_block / df_block
                ms_error <- ss_error / df_error

                var_error <- ms_error
                var_treat <- max(0, (ms_treat - ms_error) / repli)
                var_block <- max(0, (ms_block - ms_error) / genotype)

                total_var <- var_treat + var_block + var_error
                if (total_var > 0) {
                  shrinkage <- var_treat / total_var
                } else {
                  shrinkage <- 0.5
                }
              } else {
                shrinkage <- 0.5
              }

              treat_missing <- gen_idx[missing_idx]
              block_missing <- rep_idx[missing_idx]

              treat_mean_vals <- sapply(treat_missing, function(g) {
                if (g <= length(treat_means) && g > 0) treat_means[g] else grand_mean
              })
              block_mean_vals <- sapply(block_missing, function(r) {
                if (r <= length(block_means) && r > 0) block_means[r] else grand_mean
              })

              predicted <- treat_mean_vals + block_mean_vals - grand_mean
              data_imputed[missing_idx, col] <- shrinkage * predicted + (1 - shrinkage) * grand_mean
            } else {
              col_complete <- col_idx[complete_idx]

              col_sums <- grouped_sums(complete_data, col_complete)[, 1]
              col_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), col_complete)[, 1]
              unique_cols <- sort(unique(col_complete))
              col_means_vec <- col_sums / pmax(col_counts, 1)
              col_means <- numeric(genotype) # LSD uses genotype for ncol_blocks
              col_means[unique_cols] <- col_means_vec

              ss_col <- sum(col_counts * (col_means_vec - grand_mean)^2)

              residuals <- y_complete - treat_means[treat_complete] -
                block_means[block_complete] - col_means[col_complete] + 2 * grand_mean
              ss_error <- sum(residuals^2)

              df_treat <- genotype - 1
              df_row <- genotype - 1
              df_col <- genotype - 1
              df_error <- (genotype - 1) * (genotype - 2)

              if (df_error > 0 && length(complete_idx) > (3 * genotype - 2)) {
                ms_treat <- ss_treat / df_treat
                ms_row <- ss_block / df_row
                ms_col <- ss_col / df_col
                ms_error <- ss_error / df_error

                var_error <- ms_error
                var_treat <- max(0, (ms_treat - ms_error) / genotype)
                var_row <- max(0, (ms_row - ms_error) / genotype)
                var_col <- max(0, (ms_col - ms_error) / genotype)

                total_var <- var_treat + var_row + var_col + var_error
                if (total_var > 0) {
                  shrinkage <- var_treat / total_var
                } else {
                  shrinkage <- 0.33
                }
              } else {
                shrinkage <- 0.33
              }

              treat_missing <- gen_idx[missing_idx]
              row_missing <- rep_idx[missing_idx]
              col_missing <- col_idx[missing_idx]

              treat_mean_vals <- sapply(treat_missing, function(g) {
                if (g <= length(treat_means) && g > 0) treat_means[g] else grand_mean
              })
              row_mean_vals <- sapply(row_missing, function(r) {
                if (r <= length(block_means) && r > 0) block_means[r] else grand_mean
              })
              col_mean_vals <- sapply(col_missing, function(c) {
                if (c <= length(col_means) && c > 0) col_means[c] else grand_mean
              })

              predicted <- treat_mean_vals + row_mean_vals + col_mean_vals - 2 * grand_mean
              data_imputed[missing_idx, col] <- shrinkage * predicted + (1 - shrinkage) * grand_mean
            }
          }

          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
  } else {



    for (col in seq_len(ncols)) {
      missing_idx <- which(!is.finite(data_imputed[, col]))
      n_missing <- length(missing_idx)

      if (n_missing > 0) {
        col_mean <- mean(data_imputed[, col], na.rm = TRUE)
        data_imputed[missing_idx, col] <- col_mean

        for (iter in seq_len(max_iter)) {
          old_values <- data_imputed[missing_idx, col]

          complete_idx <- which(is.finite(data_mat[, col]))

          if (length(complete_idx) > 3) {
            covariate_cols <- setdiff(seq_len(ncols), col)

            valid_covariates <- integer(0)
            for (cov_col in covariate_cols) {
              if (all(is.finite(data_imputed[complete_idx, cov_col]))) {
                valid_covariates <- c(valid_covariates, cov_col)
              }
            }

            if (length(valid_covariates) > 0) {
              treat_complete <- gen_idx[complete_idx]
              treat_missing <- gen_idx[missing_idx]
              block_complete <- rep_idx[complete_idx]
              block_missing <- rep_idx[missing_idx]
              y_complete <- data_imputed[complete_idx, col]

              n_complete <- length(complete_idx)

              X_treat <- matrix(0, n_complete, genotype - 1)
              for (g in seq_len(genotype - 1)) {
                X_treat[treat_complete == g, g] <- 1
              }
              X_treat[treat_complete == genotype, ] <- -1

              X_block <- matrix(0, n_complete, repli - 1)
              for (r in seq_len(repli - 1)) {
                X_block[block_complete == r, r] <- 1
              }
              X_block[block_complete == repli, ] <- -1

              if (design_type == "LSD") {
                col_complete <- col_idx[complete_idx]
                col_missing <- col_idx[missing_idx]
                X_col <- matrix(0, n_complete, ncol_blocks - 1)
                for (c in seq_len(ncol_blocks - 1)) {
                  X_col[col_complete == c, c] <- 1
                }
                X_col[col_complete == ncol_blocks, ] <- -1
              }

              X_cov <- matrix(0, n_complete, length(valid_covariates))
              cov_means <- numeric(length(valid_covariates))
              for (k in seq_along(valid_covariates)) {
                cov_col <- valid_covariates[k]
                cov_values <- data_imputed[complete_idx, cov_col]
                cov_means[k] <- mean(cov_values)
                X_cov[, k] <- cov_values - cov_means[k]
              }

              if (design_type == "LSD") {
                X <- cbind(1, X_treat, X_block, X_col, X_cov)
              } else {
                X <- cbind(1, X_treat, X_block, X_cov)
              }

              qr_fit <- qr(X)
              if (qr_fit$rank == ncol(X)) {
                beta <- qr.coef(qr_fit, y_complete)

                X_missing <- matrix(0, n_missing, ncol(X))

                X_missing[, 1] <- 1

                for (g in seq_len(genotype - 1)) {
                  X_missing[treat_missing == g, 1 + g] <- 1
                }
                X_missing[treat_missing == genotype, 2:genotype] <- -1

                block_start <- genotype
                for (r in seq_len(repli - 1)) {
                  X_missing[block_missing == r, block_start + r] <- 1
                }
                X_missing[block_missing == repli, (block_start + 1):(block_start + repli - 1)] <- -1

                if (design_type == "LSD") {
                  col_start <- genotype + repli - 1
                  for (c in seq_len(ncol_blocks - 1)) {
                    X_missing[col_missing == c, col_start + c] <- 1
                  }
                  X_missing[col_missing == ncol_blocks, (col_start + 1):(col_start + ncol_blocks - 1)] <- -1
                  cov_start <- genotype + repli + ncol_blocks - 2
                } else {
                  cov_start <- genotype + repli - 1
                }

                for (k in seq_along(valid_covariates)) {
                  cov_col <- valid_covariates[k]
                  cov_vals <- data_imputed[missing_idx, cov_col]
                  valid_mask <- is.finite(cov_vals)
                  X_missing[valid_mask, cov_start + k] <- cov_vals[valid_mask] - cov_means[k]
                }

                data_imputed[missing_idx, col] <- X_missing %*% beta
              }
            } else {
              complete_data <- matrix(data_imputed[complete_idx, col], ncol = 1)
              treat_complete <- gen_idx[complete_idx]
              block_complete <- rep_idx[complete_idx]

              treat_sums <- grouped_sums(complete_data, treat_complete)[, 1]
              treat_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), treat_complete)[, 1]
              unique_treats <- sort(unique(treat_complete))
              treat_means_vec <- treat_sums / pmax(treat_counts, 1)
              treat_means <- numeric(genotype)
              treat_means[unique_treats] <- treat_means_vec

              block_sums <- grouped_sums(complete_data, block_complete)[, 1]
              block_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), block_complete)[, 1]
              unique_blocks <- sort(unique(block_complete))
              block_means_vec <- block_sums / pmax(block_counts, 1)
              block_means <- numeric(repli)
              block_means[unique_blocks] <- block_means_vec

              grand_mean <- cpp_grand_means(complete_data)[1]

              if (design_type == "LSD") {
                col_complete <- col_idx[complete_idx]
                col_sums <- grouped_sums(complete_data, col_complete)[, 1]
                col_counts <- grouped_sums(matrix(rep(1, length(complete_idx)), ncol = 1), col_complete)[, 1]
                unique_cols <- sort(unique(col_complete))
                col_means_vec <- col_sums / pmax(col_counts, 1)
                col_means <- numeric(ncol_blocks)
                col_means[unique_cols] <- col_means_vec

                treat_missing <- gen_idx[missing_idx]
                row_missing <- rep_idx[missing_idx]
                col_missing <- col_idx[missing_idx]

                treat_eff <- sapply(treat_missing, function(g) {
                  if (g <= length(treat_means) && g > 0) treat_means[g] - grand_mean else 0
                })
                row_eff <- sapply(row_missing, function(r) {
                  if (r <= length(block_means) && r > 0) block_means[r] - grand_mean else 0
                })
                col_eff <- sapply(col_missing, function(c) {
                  if (c <= length(col_means) && c > 0) col_means[c] - grand_mean else 0
                })

                data_imputed[missing_idx, col] <- grand_mean + treat_eff + row_eff + col_eff
              } else {
                treat_missing <- gen_idx[missing_idx]
                block_missing <- rep_idx[missing_idx]

                treat_eff <- sapply(treat_missing, function(g) {
                  if (g <= length(treat_means) && g > 0) treat_means[g] - grand_mean else 0
                })
                block_eff <- sapply(block_missing, function(r) {
                  if (r <= length(block_means) && r > 0) block_means[r] - grand_mean else 0
                })

                data_imputed[missing_idx, col] <- grand_mean + treat_eff + block_eff
              }
            }
          }

          max_change <- max(abs(data_imputed[missing_idx, col] - old_values))
          if (max_change < tolerance) break
        }
      }
    }
  }

  data_imputed
}
