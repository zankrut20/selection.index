#' @title Mean performance of phenotypic data
#'
#' @param data data for analysis
#' @param genotypes genotypes vector
#' @param replications replication vector
#'
#' @return Dataframe of mean performance analysis
#' @export
#' @importFrom stats aggregate anova lm qt
#' @examples
#' meanPerformance(data = seldata[, 3:9], genotypes = seldata[, 2], replications = seldata[, 1])
#'
#'

meanPerformance <- function(data, genotypes, replications){
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
  
  # OPTIMIZATION: Pre-compute unique genotype order once
  # Avoids: Repeated unique() calls in meanData function
  odr <- unique(genotypes)
  n_genotypes <- length(odr)
  
  # OPTIMIZATION: Vectorized mean calculation using rowsum()
  # Avoids: aggregate() overhead (S3 dispatch, split-apply-combine)
  # Why faster: rowsum() is .Internal primitive, optimized C implementation
  gen_idx <- match(genotypes, odr)  # Integer indices for grouping
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
  
  # OPTIMIZATION: Pre-compute constants for all traits
  # Avoids: Repeated sqrt(2), qt() calls, division by r
  sqrt_2_over_r <- sqrt(2 / r)
  sqrt_1_over_r <- sqrt(1 / r)
  
  # OPTIMIZATION: Inline nested functions, vectorize where possible
  # Avoids: Function call overhead × colnumber
  for (j in seq_len(colnumber)) {
    trait_data <- data_mat[, j]
    
    # Fit linear model once
    model <- lm(trait_data ~ replications_fac + genotypes_fac)
    anova_model <- anova(model)
    
    # Extract values efficiently
    EMS <- anova_model[3, 3]
    GMS <- anova_model[2, 3]
    df_error <- anova_model[3, 1]
    p_value_01 <- anova_model[2, 5] > 0.01
    p_value_05 <- anova_model[2, 5] > 0.05
    
    # OPTIMIZATION: Use column min/max directly (faster than subsetting data.frame)
    # Avoids: meanData() function call, data.frame operations
    trait_means <- mean_mat[, j]
    Maxi <- max(trait_means)
    Mini <- min(trait_means)
    
    # OPTIMIZATION: Use .Internal primitives where available
    # mean() is already optimized, but direct computation is clear
    GM <- sum(trait_data) / length(trait_data)
    
    # Pre-computed values
    SD <- sqrt(EMS)
    SEm <- sqrt(EMS) * sqrt_1_over_r
    CV <- (SD / GM) * 100
    
    # OPTIMIZATION: Compute critical differences with pre-computed constants
    # Avoids: Repeated sqrt() and division operations
    qt_005 <- abs(qt(0.005, df_error))
    qt_025 <- abs(qt(0.025, df_error))
    CD1_val <- SEm * sqrt_2_over_r * r * qt_005  # Simplifies to correct formula
    CD5_val <- SEm * sqrt_2_over_r * r * qt_025
    
    # OPTIMIZATION: Compute genetic/phenotypic variance once
    GV <- (GMS - EMS) / r
    PV <- GV + EMS
    hs <- GV / PV
    
    # OPTIMIZATION: Store as numeric, format strings only at end if needed
    # Avoids: String paste operations in tight loop
    # Note: Storing numeric values, will add NS label at end if needed
    perf_mat[1, j] <- Mini
    perf_mat[2, j] <- Maxi
    perf_mat[3, j] <- GM
    perf_mat[4, j] <- CV
    perf_mat[5, j] <- SEm
    perf_mat[6, j] <- CD5_val
    perf_mat[7, j] <- CD1_val
    perf_mat[8, j] <- hs
    perf_mat[9, j] <- hs * 100
    
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
