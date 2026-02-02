#' Genotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments
#' @param replication vector containing replication
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return A Genotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' gen.varcov(data=seldata[,3:9], genotypes=seldata$treat,replication=seldata$rep)
gen.varcov<- function (data, genotypes, replication, method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"))
{
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
  
  # OPTIMIZATION: Pre-allocate result matrix (not growing vector)
  # Avoids: Memory reallocation on every c(x, new_value) - O(n²) copies
  # Why faster: Single allocation, direct indexing, no memory churn
  genetic.cov <- matrix(0, nrow = colnumber, ncol = colnumber,
                       dimnames = list(headings, headings))
  
  # DESIGN ENGINE: Use modular RCBD calculations
  # Eliminates ~60 lines of repeated design calculation code
  # Centralizes correction factor, sums of products, mean products computations
  for (i in seq_len(colnumber)) {
    trait1 <- data_mat[, i]
    
    for (j in seq_len(colnumber)) {
      trait2 <- data_mat[, j]
      
      # Single call to design engine replaces ~15 lines of manual calculations
      design_stats <- rcbd.design(trait1, trait2, gen_idx, rep_idx, 
                                   calc_type = "mean_products")
      
      # Genotypic covariance = (GMP - EMP) / r
      genetic.cov[i, j] <- (design_stats$GMP - design_stats$EMP) / design_stats$n_replications
    }
  }
  
  return(genetic.cov)
}
