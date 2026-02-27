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

  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }

  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"

  colnumber <- ncol(data_mat)
  col_names <- colnames(data)

  genotypes_fac <- factor(genotypes, levels = unique(genotypes))
  replications_fac <- as.factor(replications)
  r <- nlevels(replications_fac)

  if (design_type == "LSD") {
    columns_fac <- as.factor(columns)
  }

  if (design_type == "SPD") {
    main_plots_fac <- as.factor(main_plots)
  }

  odr <- levels(genotypes_fac)


  if (any(!is.finite(data_mat))) {
    method_provided <- !missing(method)
    method <- match.arg(method)

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

  gen_idx <- as.integer(genotypes_fac)
  rep_idx <- as.integer(replications_fac)
  col_idx <- if (design_type == "LSD") as.integer(columns_fac) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots_fac) else NULL

  mean_mat <- rowsum(data_mat, gen_idx, reorder = FALSE) / tabulate(gen_idx)
  mean_mat <- round(mean_mat, 4)
  meandf <- data.frame(
    Genotypes = odr, mean_mat,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  colnames(meandf) <- c("Genotypes", col_names)

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

  GMS_vec <- anova_result$GMS
  EMS_vec <- anova_result$EMS

  df_genotype <- anova_result$DFG
  df_error <- anova_result$DFE

  n_main <- anova_result$n_main

  perf_mat <- matrix(0, nrow = 9, ncol = colnumber)
  perf_labels <- matrix("", nrow = 9, ncol = colnumber) # For storing NS labels per trait

  for (j in seq_len(colnumber)) {
    trait_data <- data_mat[, j]

    GMS <- GMS_vec[j]
    EMS <- EMS_vec[j]



    F_stat <- if (!is.na(GMS) && !is.na(EMS) && EMS > 0) GMS / EMS else NA_real_
    p_value <- if (!is.na(F_stat) && !is.na(df_genotype) && !is.na(df_error)) {
      pf(F_stat, df_genotype, df_error, lower.tail = FALSE)
    } else {
      NA_real_
    }

    p_value_01 <- if (!is.na(p_value)) p_value > 0.01 else FALSE
    p_value_05 <- if (!is.na(p_value)) p_value > 0.05 else FALSE

    trait_means <- mean_mat[, j]
    Maxi <- max(trait_means, na.rm = TRUE)
    Mini <- min(trait_means, na.rm = TRUE)

    GM <- mean(trait_data, na.rm = TRUE)

    if (design_type == "SPD") {
      SEm <- if (!is.na(EMS) && r > 0 && n_main > 0) sqrt(EMS / (r * n_main)) else NA_real_

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
      SEm <- if (!is.na(EMS) && r > 0) sqrt(EMS / r) else NA_real_

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

    CV <- if (!is.na(EMS) && !is.na(GM) && GM != 0) {
      (sqrt(EMS) / GM) * 100
    } else {
      NA_real_
    }

    if (design_type == "SPD") {
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

    PV <- if (!is.na(GV) && !is.na(EMS)) GV + EMS else NA_real_

    hs <- if (!is.na(GV) && !is.na(PV) && PV > 0) GV / PV else NA_real_

    perf_mat[1, j] <- Mini
    perf_mat[2, j] <- Maxi
    perf_mat[3, j] <- GM
    perf_mat[4, j] <- CV
    perf_mat[5, j] <- SEm
    perf_mat[6, j] <- CD5_val
    perf_mat[7, j] <- CD1_val
    perf_mat[8, j] <- hs
    perf_mat[9, j] <- if (!is.na(hs)) hs * 100 else NA_real_

    if (p_value_05) {
      perf_labels[6, j] <- " NS"
    }
    if (p_value_01) {
      perf_labels[7, j] <- " NS"
    }
  }

  perf_mat <- round(perf_mat, 4)

  perf_df <- as.data.frame(perf_mat, stringsAsFactors = FALSE)
  colnames(perf_df) <- col_names

  for (j in seq_len(colnumber)) {
    if (nzchar(perf_labels[6, j])) {
      perf_df[6, j] <- paste0(perf_df[6, j], perf_labels[6, j])
    }
    if (nzchar(perf_labels[7, j])) {
      perf_df[7, j] <- paste0(perf_df[7, j], perf_labels[7, j])
    }
  }

  row_labels <- c(
    "Min", "Max", "GM", "CV (%)", "SEm", "CD 5%", "CD 1%",
    "Heritability", "Heritability(%)"
  )
  perf_df <- data.frame(
    Genotypes = row_labels, perf_df,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  result <- rbind(meandf, perf_df)

  rownames(result) <- NULL

  result
}
