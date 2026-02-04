#' @title Mean performance of phenotypic data
#'
#' @param data data for analysis
#' @param genotypes genotypes vector
#' @param replications replication vector
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param design_type experimental design type: "RCBD" (default) or "LSD" (Latin Square)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return Dataframe of mean performance analysis
#' @export
#' @importFrom stats qt pf
#' @examples
#' meanPerformance(data = seldata[, 3:9], genotypes = seldata[, 2], replications = seldata[, 1])
#'
#'

meanPerformance <- function(data, genotypes, replications, columns = NULL, 
                           design_type = c("RCBD", "LSD"),
                           method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett")){
  design_type <- match.arg(design_type)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  
  # OPTIMIZATION: Convert to numeric matrix once (avoid repeated conversions)
  # Avoids: sapply(data, as.numeric) + as.list() overhead in nested functions
  # Why faster: Single type coercion, direct matrix column access
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  
  colnumber <- ncol(data_mat)
  col_names <- colnames(data)
  
  # OPTIMIZATION: Factor conversion outside loops
  # Avoids: Repeated as.factor() calls (was called colnumber times)
  # Why faster: Factor creation is expensive (level sorting, attribute setup)
  genotypes_fac <- as.factor(genotypes)
  replications_fac <- as.factor(replications)
  r <- nlevels(replications_fac)
  
  if (design_type == "LSD") {
    columns_fac <- as.factor(columns)
  }
  
  # OPTIMIZATION: Pre-compute unique genotype order once
  # Avoids: Repeated unique() calls in meanData function
  odr <- unique(genotypes)
  n_genotypes <- length(odr)
  
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
              call. = FALSE)
    }
    
    gen_idx <- as.integer(genotypes_fac)
    rep_idx <- as.integer(replications_fac)
    col_idx <- if (design_type == "LSD") as.integer(columns_fac) else NULL
    data_mat <- missingValueEstimation(data_mat, gen_idx, rep_idx, col_idx, design_type, method)
  }
  
  # Convert to integer indices for design.stats engine
  gen_idx <- as.integer(genotypes_fac)
  rep_idx <- as.integer(replications_fac)
  col_idx <- if (design_type == "LSD") as.integer(columns_fac) else NULL
  
  # OPTIMIZATION: Vectorized mean calculation using rowsum()
  # Avoids: aggregate() overhead (S3 dispatch, split-apply-combine)
  # Why faster: rowsum() is .Internal primitive, optimized C implementation
  mean_mat <- rowsum(data_mat, gen_idx, reorder = FALSE) / tabulate(gen_idx)
  mean_mat <- round(mean_mat, 4)
  meandf <- data.frame(Genotypes = odr, mean_mat, 
                       stringsAsFactors = FALSE, check.names = FALSE)
  colnames(meandf) <- c("Genotypes", col_names)
  
  # OPTIMIZATION: Pre-allocate performance matrix (not growing list)
  # Avoids: Memory reallocation in list growth
  # Why faster: Single allocation, direct column assignment
  perf_mat <- matrix(0, nrow = 9, ncol = colnumber)
  perf_labels <- matrix(character(9), ncol = 1)  # For storing NS labels
  
  # OPTIMIZATION: Use design.stats engine for ANOVA calculations
  # Avoids: Repeated lm/anova overhead, centralizes design calculations
  # Supports both RCBD and LSD designs
  for (j in seq_len(colnumber)) {
    trait_data <- data_mat[, j]
    
    # Use design.stats engine to get mean squares and degrees of freedom
    if (design_type == "RCBD") {
      design_result <- design.stats(trait_data, trait_data, gen_idx, rep_idx,
                                   design_type = "RCBD", calc_type = "all")
    } else {
      design_result <- design.stats(trait_data, trait_data, gen_idx, rep_idx, col_idx,
                                   design_type = "LSD", calc_type = "all")
    }
    
    # Extract design statistics
    EMS <- design_result$EMP  # Error Mean Product (for variance, this is MSE)
    GMS <- design_result$GMP  # Genotype Mean Product (for variance, this is MSG)
    df_error <- design_result$DFE
    
    # For significance testing, we need F-statistic and p-value
    # F = GMS / EMS
    F_stat <- if (!is.na(GMS) && !is.na(EMS) && EMS > 0) GMS / EMS else NA_real_
    p_value <- if (!is.na(F_stat) && !is.na(design_result$DFG) && !is.na(df_error)) {
      pf(F_stat, design_result$DFG, df_error, lower.tail = FALSE)
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
    
    # Coefficient of Variation: CV% = (sqrt(MSE) / grand_mean) * 100
    CV <- if (!is.na(EMS) && !is.na(GM) && GM != 0) {
      (sqrt(EMS) / GM) * 100
    } else {
      NA_real_
    }
    
    # Genetic Variance: GV = (MSG - MSE) / r (guard against negative)
    GV <- if (!is.na(GMS) && !is.na(EMS) && r > 0) {
      max(0, (GMS - EMS) / r)
    } else {
      NA_real_
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
    
    # Store significance labels separately (avoid mixing types)
    if (p_value_05) {
      perf_labels[6] <- "NS"
    }
    if (p_value_01) {
      perf_labels[7] <- "NS"
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
  
  # Add NS labels where needed (paste only if non-empty)
  for (i in c(6, 7)) {
    if (nzchar(perf_labels[i])) {
      perf_df[i, ] <- vapply(perf_df[i, ], function(x) {
        paste(x, perf_labels[i])
      }, character(1), USE.NAMES = FALSE)
    }
  }
  
  # Add row names as first column
  row_labels <- c("Min", "Max", "GM", "CV (%)", "SEm", "CD 5%", "CD 1%", 
                  "Heritability", "Heritability(%)")
  perf_df <- data.frame(Genotypes = row_labels, perf_df, 
                        stringsAsFactors = FALSE, check.names = FALSE)
  
  # OPTIMIZATION: Single rbind at end (not incremental)
  # Avoids: Growing data.frame row-by-row
  # Why faster: Pre-allocated result, single memory operation
  result <- rbind(meandf, perf_df)
  
  return(result)
}
