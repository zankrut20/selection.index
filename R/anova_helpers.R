#' Calculate ANOVA Components (Internal Helper - DEPRECATED)
#'
#' @description
#' **DEPRECATED**: This function is deprecated and should not be used.
#' Use \code{\link{design_stats_api}} instead, which provides centralized
#' ANOVA computation through the design_stats.R engine.
#'
#' Internal function to compute ANOVA mean squares using math primitives.
#' Kept for backward compatibility only.
#'
#' @param data_mat Numeric matrix of trait data
#' @param gen_idx Integer vector of genotype indices
#' @param rep_idx Integer vector of replication indices
#' @param col_idx Integer vector of column indices (for LSD, optional)
#' @param main_idx Integer vector of main plot indices (for SPD, optional)
#' @param design_type Integer design code: 1=RCBD, 2=LSD, 3=SPD
#'
#' @return List with components:
#'   \item{GMS}{Genotype mean squares vector}
#'   \item{EMS}{Error mean squares vector}
#'   \item{MSG}{Genotype mean squares matrix (for variance-covariance)}
#'   \item{MSE}{Error mean squares matrix (for variance-covariance)}
#'   \item{df_gen}{Genotype degrees of freedom}
#'   \item{df_error}{Error degrees of freedom}
#'
#' @keywords internal
#' @noRd
.calculate_anova <- function(data_mat, gen_idx, rep_idx,
                             col_idx = NULL, main_idx = NULL,
                             design_type = 1L) {
  # DEPRECATION WARNING
  warning(
    ".calculate_anova() is deprecated. Use design_stats_api() instead.\n",
    "design_stats.R is now the single engine for ANOVA computations.",
    call. = FALSE
  )

  n_traits <- ncol(data_mat)
  n_obs <- nrow(data_mat)

  # Grouped statistics using math primitive
  gen_sums <- grouped_sums(data_mat, gen_idx)
  rep_sums <- grouped_sums(data_mat, rep_idx)

  # Calculate counts for each group
  gen_counts <- as.integer(table(gen_idx))
  rep_counts <- as.integer(table(rep_idx))

  # Correction factor using math primitive
  total_sums <- colSums(data_mat)
  CF <- correction_factor(total_sums, n_obs)

  # Sum of products using math primitives
  TSP <- total_sum_of_products(data_mat, CF)
  GSP <- grouped_sum_of_products(gen_sums, gen_counts, CF)
  RSP <- grouped_sum_of_products(rep_sums, rep_counts, CF)

  # Design-specific calculations
  if (design_type == 1L) { # RCBD
    n_gen <- nrow(gen_sums)
    n_rep <- nrow(rep_sums)
    df_gen <- n_gen - 1
    df_rep <- n_rep - 1
    df_error <- df_gen * df_rep

    ESP <- TSP - GSP - RSP
    MSG <- mean_squares(GSP, df_gen)
    MSR <- mean_squares(RSP, df_rep)
    MSE <- mean_squares(ESP, df_error)

    # For SPD compatibility
    EMS_MAIN_vec <- rep(NA_real_, n_traits)
    DFE_MAIN <- NA_integer_
    n_main <- NA_integer_
  } else if (design_type == 2L) { # LSD
    col_sums <- grouped_sums(data_mat, col_idx)
    col_counts <- as.integer(table(col_idx))

    n_gen <- nrow(gen_sums)
    n_rep <- nrow(rep_sums)
    n_col <- nrow(col_sums)

    df_gen <- n_gen - 1
    df_rep <- n_rep - 1
    df_col <- n_col - 1
    df_error <- (n_gen - 1) * (n_rep - 1) - (n_col - 1)

    CSP <- grouped_sum_of_products(col_sums, col_counts, CF)
    ESP <- TSP - GSP - RSP - CSP

    MSG <- mean_squares(GSP, df_gen)
    MSR <- mean_squares(RSP, df_rep)
    MSC <- mean_squares(CSP, df_col)
    MSE <- mean_squares(ESP, df_error)

    # For SPD compatibility
    EMS_MAIN_vec <- rep(NA_real_, n_traits)
    DFE_MAIN <- NA_integer_
    n_main <- NA_integer_
  } else if (design_type == 3L) { # SPD
    main_sums <- grouped_sums(data_mat, main_idx)
    main_counts <- as.integer(table(main_idx))

    n_gen <- nrow(gen_sums)
    n_rep <- nrow(rep_sums)
    n_main <- nrow(main_sums)

    df_gen <- n_gen - 1
    df_main <- n_main - 1
    df_rep <- n_rep - 1
    df_error1 <- (n_gen - 1) * (n_main - 1)
    df_error2 <- n_gen * (n_rep - 1)

    MSP <- grouped_sum_of_products(main_sums, main_counts, CF)
    ESP1 <- GSP - MSP
    ESP2 <- TSP - GSP - RSP

    MSG <- mean_squares(GSP, df_gen)
    MSM <- mean_squares(MSP, df_main)
    MSE1 <- mean_squares(ESP1, df_error1)
    MSE2 <- mean_squares(ESP2, df_error2)

    # For SPD, use MSE2 as the primary error term
    MSE <- MSE2
    df_error <- df_error2

    # Extract EMS_MAIN for SPD
    EMS_MAIN_vec <- diag(MSE1)
    DFE_MAIN <- df_error1
  }

  # Extract diagonal elements for mean_performance compatibility
  GMS_vec <- diag(MSG)
  EMS_vec <- diag(MSE)

  # Return results in same format as cpp_anova_iterator
  list(
    GMS = GMS_vec, # Vector for mean_performance
    EMS = EMS_vec, # Vector for mean_performance
    EMS_MAIN = EMS_MAIN_vec, # Vector for SPD
    DFG = df_gen, # Degrees of freedom genotype
    DFE = df_error, # Degrees of freedom error
    DFE_MAIN = DFE_MAIN, # Degrees of freedom error main (SPD)
    n_rep = length(unique(rep_idx)),
    n_gen = length(unique(gen_idx)),
    n_main = ifelse(design_type == 3L, n_main, NA_integer_),
    MSG = MSG, # Matrix for variance-covariance
    MSE = MSE # Matrix for variance-covariance
  )
}
