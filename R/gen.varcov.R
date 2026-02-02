#' Genotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments
#' @param replication vector containing replication/blocks (RCBD) or rows (LSD)
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param design_type experimental design type: "RCBD" (default) or "LSD" (Latin Square)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return A Genotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' # RCBD example
#' gen.varcov(data=seldata[,3:9], genotypes=seldata$treat, replication=seldata$rep)
#' 
#' # Latin Square Design example (requires columns parameter)
#' # gen.varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat, 
#' #           replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")
gen.varcov<- function (data, genotypes, replication, columns = NULL, design_type = c("RCBD", "LSD"), 
                       method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"))
{
  design_type <- match.arg(design_type)
  
  # OPTIMIZATION: Single matrix conversion with storage.mode assignment
  # Avoids: (1) data.frame->list conversion overhead, (2) repeated as.numeric() per column
  # Why faster: Direct storage type coercion in C, no intermediate structures
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  
  colnumber <- ncol(data_mat)
  headings <- colnames(data)
  
  # OPTIMIZATION: Convert factors once outside loops (not colnumber² times)
  # Avoids: Redundant as.factor() calls and nlevels() computations
  # Why faster: Factor conversion is expensive (level sorting, attribute creation)
  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)
  repli <- nlevels(replication)
  genotype <- nlevels(genotypes)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  if (design_type == "LSD") {
    columns <- as.factor(columns)
  }
  
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
    
    gen_idx <- as.integer(genotypes)
    rep_idx <- as.integer(replication)
    data_mat <- missingValueEstimation(data_mat, gen_idx, rep_idx, method)
  }
  
  # OPTIMIZATION: Convert factors to integer indices for design engine
  # Avoids: Factor level lookups in internal calculations
  # Why faster: Integer indexing is primitive operation
  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
  
  # OPTIMIZATION: Pre-allocate result matrix (not growing vector)
  # Avoids: Memory reallocation on every c(x, new_value) - O(n²) copies
  # Why faster: Single allocation, direct indexing, no memory churn
  genetic.cov <- matrix(0, nrow = colnumber, ncol = colnumber,
                       dimnames = list(headings, headings))
  
  # DESIGN ENGINE: Use modular design calculations (RCBD or LSD)
  # Eliminates ~60 lines of repeated design calculation code
  # Centralizes correction factor, sums of products, mean products computations
  for (i in seq_len(colnumber)) {
    trait1 <- data_mat[, i]
    
    for (j in seq_len(colnumber)) {
      trait2 <- data_mat[, j]
      
      # Single call to design engine replaces ~15 lines of manual calculations
      if (design_type == "RCBD") {
        design_stats <- rcbd.design(trait1, trait2, gen_idx, rep_idx, 
                                     design_type = "RCBD", calc_type = "mean_products")
        # Genotypic covariance = (GMP - EMP) / r
        genetic.cov[i, j] <- (design_stats$GMP - design_stats$EMP) / design_stats$n_replications
      } else {
        design_stats <- rcbd.design(trait1, trait2, gen_idx, rep_idx, col_idx,
                                     design_type = "LSD", calc_type = "mean_products")
        # For LSD: Genotypic covariance = (GMP - EMP) / t
        genetic.cov[i, j] <- (design_stats$GMP - design_stats$EMP) / design_stats$n_genotypes
      }
    }
  }
  
  return(genetic.cov)
}
