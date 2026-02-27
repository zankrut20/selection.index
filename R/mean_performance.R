#' @title Mean performance of phenotypic data
#'
#' @param data data for analysis
#' @param genotypes genotypes vector (sub-plot treatments in SPD)
#' @param replications replication vector
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param main_plots vector containing main plot treatments (required for Split Plot Design only)
#' @param design_type experimental design type: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return Dataframe of mean performance analysis
#' @export
#' @importFrom stats qt pf
#' @examples
#' mean_performance(data = seldata[, 3:9], genotypes = seldata[, 2], replications = seldata[, 1])
#'
mean_performance <- function(data, genotypes, replications, columns = NULL, main_plots = NULL,
                             design_type = c("RCBD", "LSD", "SPD"),
                             method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")) {
  design_type <- match.arg(design_type)

  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }

  # Validate Split Plot Design requirements
  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }

  # OPTIMIZATION: Convert to numeric matrix once (avoid repeated conversions)
  # Avoids: sapply(data, as.numeric) + as.list() overhead in nested functions
  # Why faster: Single type coercion, direct matrix column access
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"

  colnumber <- ncol(data_mat)
  col_names <- colnames(data)

  # OPTIMIZATION: Factor conversion outside loops (ensures genotypes/replications are categorical)
  # Avoids: Repeated as.factor() calls (was called colnumber times)
  # Why faster: Factor creation is expensive (level sorting, attribute setup)
  # Ensure genotypes and replications are always factors, even if input is numeric
  # Use factor() with ordered unique levels to preserve input order (G1, G2, G3 not G1, G10, G11)
  genotypes_fac <- factor(genotypes, levels = unique(genotypes))
  replications_fac <- as.factor(replications)
  r <- nlevels(replications_fac)

  if (design_type == "LSD") {
    columns_fac <- as.factor(columns)
  }

  if (design_type == "SPD") {
    main_plots_fac <- as.factor(main_plots)
  }

  # OPTIMIZATION: Pre-compute unique genotype order once from factor levels
  # Avoids: Repeated unique() calls in meanData function
  # Uses factor levels to preserve categorical nature of genotypes
  odr <- levels(genotypes_fac)


  # MISSING VALUE HANDLING: Use modular engine for imputation
  # Only process method parameter if missing values are detected
  if (any(!is.finite(data_mat))) {
    # Check if user explicitly provided a method
    method_provided <- !missing(method)
    method <- match.arg(method)

    # Warn user if they have missing values but didn't explicitly specify a method
    if (!method_provided) {
      warning("Missing values detected in data. Using default method 'REML' for imputation. ",
        "Consider explicitly specifying method: 'REML', 'Yates', 'Healy', 'Regression', 'Mean', or 'Bartlett'.",
        call. = FALSE
      )
    }

    gen_idx <- as.integer(genotypes_fac)
    rep_idx <- as.integer(replications_fac)
    col_idx <- if (design_type == "LSD") as.integer(columns_fac) else NULL
    main_idx <- if (design_type == "SPD") as.integer(main_plots_fac) else NULL
    data_mat <- missing_value_estimation(data_mat, gen_idx, rep_idx, col_idx, main_idx, design_type, method)
  }

  # Convert to integer indices for design_stats engine
  gen_idx <- as.integer(genotypes_fac)
  rep_idx <- as.integer(replications_fac)
  col_idx <- if (design_type == "LSD") as.integer(columns_fac) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots_fac) else NULL

  # OPTIMIZATION: Vectorized mean calculation using rowsum()
  # Avoids: aggregate() overhead (S3 dispatch, split-apply-combine)
  # Why faster: rowsum() is .Internal primitive, optimized C implementation
  mean_mat <- rowsum(data_mat, gen_idx, reorder = FALSE) / tabulate(gen_idx)
  mean_mat <- round(mean_mat, 4)
  meandf <- data.frame(
    Genotypes = odr, mean_mat,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(meandf) <- c("Genotypes", col_names)

  # CENTRALIZED: Use design_stats_api as single engine for ANOVA
  # Replaces ad-hoc vectorized ANOVA computation
  # design_stats.R is now the single source of truth for:
  # - Correction factors, sums of products, DF, MSG/MSE
  # Uses optimized design_stats() engine for all trait pairs
  design_code <- switch(design_type,
    "RCBD" = 1L,
    "LSD" = 2L,
    "SPD" = 3L
  )

  anova_result <- design_stats_api(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design_code
  )

  # Extract vectors for all traits at once
  GMS_vec <- anova_result$GMS
  EMS_vec <- anova_result$EMS

  df_genotype <- anova_result$DFG
  df_error <- anova_result$DFE

  n_main <- anova_result$n_main

  # OPTIMIZATION: Pre-allocate performance matrix (not growing list)
  # Avoids: Memory reallocation in list growth
  # Why faster: Single allocation, direct column assignment
  perf_mat <- matrix(0, nrow = 9, ncol = colnumber)
  perf_labels <- matrix("", nrow = 9, ncol = colnumber) # For storing NS labels per trait

  # Process each trait with pre-computed ANOVA statistics
  for (j in seq_len(colnumber)) {
    trait_data <- data_mat[, j]

    # Extract pre-computed statistics for this trait
    GMS <- GMS_vec[j]
    EMS <- EMS_vec[j]


    # For significance testing, we need F-statistic and p-value

    F_stat <- if (!is.na(GMS) && !is.na(EMS) && EMS > 0) GMS / EMS else NA_real_
    p_value <- if (!is.na(F_stat) && !is.na(df_genotype) && !is.na(df_error)) {
      pf(F_stat, df_genotype, df_error, lower.tail = FALSE)
    } else {
      NA_real_
    }

    p_value_01 <- if (!is.na(p_value)) p_value > 0.01 else FALSE
    p_value_05 <- if (!is.na(p_value)) p_value > 0.05 else FALSE

    # Get genotype means min/max
    trait_means <- mean_mat[, j]
    Maxi <- max(trait_means, na.rm = TRUE)
    Mini <- min(trait_means, na.rm = TRUE)

    # Grand mean
    GM <- mean(trait_data, na.rm = TRUE)

    # --- CORRECTED STATISTICAL FORMULAS ---
    # For SPD: SEm, CD use sub-plot error; accounting for nested structure
    # SEm = sqrt(MSE / (r * a)) for SPD, sqrt(MSE / r) for RCBD/LSD
    if (design_type == "SPD") {
      SEm <- if (!is.na(EMS) && r > 0 && n_main > 0) sqrt(EMS / (r * n_main)) else NA_real_

      # Critical Difference for SPD: CD = t_crit * sqrt(2 * MSE / (r * a))
      CD5_val <- if (!is.na(df_error) && !is.na(EMS) && r > 0 && n_main > 0) {
        qt(0.975, df_error) * sqrt(2 * EMS / (r * n_main))
      } else {
        NA_real_
      }

      CD1_val <- if (!is.na(df_error) && !is.na(EMS) && r > 0 && n_main > 0) {
        qt(0.995, df_error) * sqrt(2 * EMS / (r * n_main))
      } else {
        NA_real_
      }
    } else {
      # Standard Error of Mean: SEm = sqrt(MSE / r)
      SEm <- if (!is.na(EMS) && r > 0) sqrt(EMS / r) else NA_real_

      # Critical Difference (two-tailed): CD = t_crit * sqrt(2 * MSE / r)
      CD5_val <- if (!is.na(df_error) && !is.na(EMS) && r > 0) {
        qt(0.975, df_error) * sqrt(2 * EMS / r)
      } else {
        NA_real_
      }

      CD1_val <- if (!is.na(df_error) && !is.na(EMS) && r > 0) {
        qt(0.995, df_error) * sqrt(2 * EMS / r)
      } else {
        NA_real_
      }
    }

    # Coefficient of Variation: CV% = (sqrt(MSE) / grand_mean) * 100
    CV <- if (!is.na(EMS) && !is.na(GM) && GM != 0) {
      (sqrt(EMS) / GM) * 100
    } else {
      NA_real_
    }

    # Genetic Variance: For SPD: GV = (MSG - MSE) / (r * a), else (MSG - MSE) / r
    if (design_type == "SPD") {
      # n_main already extracted from anova_result above
      GV <- if (!is.na(GMS) && !is.na(EMS) && r > 0 && n_main > 0) {
        max(0, (GMS - EMS) / (r * n_main))
      } else {
        NA_real_
      }
    } else {
      GV <- if (!is.na(GMS) && !is.na(EMS) && r > 0) {
        max(0, (GMS - EMS) / r)
      } else {
        NA_real_
      }
    }

    # Phenotypic Variance: PV = GV + MSE
    PV <- if (!is.na(GV) && !is.na(EMS)) GV + EMS else NA_real_

    # Broad-sense Heritability: H² = GV / PV
    hs <- if (!is.na(GV) && !is.na(PV) && PV > 0) GV / PV else NA_real_

    # Store numeric values
    perf_mat[1, j] <- Mini
    perf_mat[2, j] <- Maxi
    perf_mat[3, j] <- GM
    perf_mat[4, j] <- CV
    perf_mat[5, j] <- SEm
    perf_mat[6, j] <- CD5_val
    perf_mat[7, j] <- CD1_val
    perf_mat[8, j] <- hs
    perf_mat[9, j] <- if (!is.na(hs)) hs * 100 else NA_real_

    # Store significance labels per trait (avoid mixing types)
    if (p_value_05) {
      perf_labels[6, j] <- " NS"
    }
    if (p_value_01) {
      perf_labels[7, j] <- " NS"
    }
  }

  # OPTIMIZATION: Round once at end (not in loop)
  # Avoids: Repeated round() function calls
  perf_mat <- round(perf_mat, 4)

  # OPTIMIZATION: Efficient data.frame construction (single call)
  # Avoids: Multiple cbind operations
  # Note: Handling NS labels - append to numeric values if needed
  perf_df <- as.data.frame(perf_mat, stringsAsFactors = FALSE)
  colnames(perf_df) <- col_names

  # Add NS labels where needed - per trait (column)
  for (j in seq_len(colnumber)) {
    if (nzchar(perf_labels[6, j])) {
      perf_df[6, j] <- paste0(perf_df[6, j], perf_labels[6, j])
    }
    if (nzchar(perf_labels[7, j])) {
      perf_df[7, j] <- paste0(perf_df[7, j], perf_labels[7, j])
    }
  }

  # Add row names as first column
  row_labels <- c(
    "Min", "Max", "GM", "CV (%)", "SEm", "CD 5%", "CD 1%",
    "Heritability", "Heritability(%)"
  )
  perf_df <- data.frame(
    Genotypes = row_labels, perf_df,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  # OPTIMIZATION: Single rbind at end (not incremental)
  # Avoids: Growing data.frame row-by-row
  # Why faster: Pre-allocated result, single memory operation
  result <- rbind(meandf, perf_df)

  # Reset row names to sequential numbers for clean output
  rownames(result) <- NULL

  result
}
