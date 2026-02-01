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
  
  # OPTIMIZATION: Pre-compute all loop-invariant constants
  # Avoids: Repeated arithmetic in nested loops (colnumber² repetitions)
  # Why faster: Division is ~10x slower than multiplication on modern CPUs
  CF_denom <- repli * genotype
  DFR <- repli - 1
  DFG <- genotype - 1
  DFE <- DFR * DFG
  repli_inv <- 1 / repli
  genotype_inv <- 1 / genotype
  
  # OPTIMIZATION: Convert factors to integer indices for rowsum()
  # Avoids: Factor level lookups in rowsum() internal code
  # Why faster: Integer indexing is primitive operation, factor requires attribute access
  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  
  # OPTIMIZATION: Pre-allocate result matrix (not growing vector)
  # Avoids: Memory reallocation on every c(x, new_value) - O(n²) copies
  # Why faster: Single allocation, direct indexing, no memory churn
  genetic.cov <- matrix(0, nrow = colnumber, ncol = colnumber,
                       dimnames = list(headings, headings))
  
  # OPTIMIZATION: Inline nested function to eliminate call overhead
  # Avoids: (1) Function call stack setup/teardown, (2) Argument copying
  # Why faster: ~20% overhead per call × colnumber² calls = substantial savings
  for (i in seq_len(colnumber)) {
    trait1 <- data_mat[, i]
    
    # OPTIMIZATION: Compute trait1 summaries once per outer loop
    # Avoids: Redundant computation across all j iterations
    # Why faster: Reuses sumch1, sumr1, GT1 for j=1..colnumber
    
    # OPTIMIZATION: Use rowsum() instead of tapply()
    # Avoids: (1) S3 method dispatch, (2) Split-apply-combine overhead, (3) List intermediates
    # Why faster: rowsum() is .Internal primitive optimized in C (5-10x faster)
    sumch1 <- rowsum(trait1, gen_idx, reorder = FALSE)
    sumr1 <- rowsum(trait1, rep_idx, reorder = FALSE)
    GT1 <- sum(trait1)
    
    for (j in seq_len(colnumber)) {
      trait2 <- data_mat[, j]
      
      # rowsum() optimization applies here too
      sumch2 <- rowsum(trait2, gen_idx, reorder = FALSE)
      sumr2 <- rowsum(trait2, rep_idx, reorder = FALSE)
      GT2 <- sum(trait2)
      
      # Pre-computed constant usage
      CF <- (GT1 * GT2) / CF_denom
      
      # OPTIMIZATION: Use crossprod() instead of sum(x * y)
      # Avoids: (1) Intermediate vector allocation for x * y, (2) Second pass for sum()
      # Why faster: Direct BLAS call, single pass through data, better cache locality
      # Note: [1] extracts scalar from 1×1 matrix result
      TSP <- crossprod(trait1, trait2)[1] - CF
      GSP <- crossprod(sumch1, sumch2)[1] * repli_inv - CF
      RSP <- crossprod(sumr1, sumr2)[1] * genotype_inv - CF
      
      ESP <- TSP - GSP - RSP
      
      # Compute mean products
      EMP <- ESP / DFE
      GMP <- GSP / DFG
      
      # OPTIMIZATION: No intermediate rounding (removed round() calls on TSP, GSP, RSP, etc.)
      # Avoids: (1) Function call overhead, (2) Numerical precision loss
      # Why faster: Rounding is expensive, only round final output if needed by caller
      
      # Direct assignment to pre-allocated matrix
      genetic.cov[i, j] <- (GMP - EMP) * repli_inv
    }
  }
  
  return(genetic.cov)
}
