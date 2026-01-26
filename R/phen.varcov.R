#' Phenotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments
#' @param replication vector containing replication
#'
#' @return A Phenotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' phen.varcov(data=seldata[,3:9], genotypes=seldata$treat,replication=seldata$rep)
phen.varcov<- function (data, genotypes, replication)
{
  # Convert to numeric matrix once - avoids repeated data.frame/list conversions
  # storage.mode assignment is faster than as.numeric() on each column
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  
  colnumber <- ncol(data_mat)
  headings <- colnames(data)
  
  # Convert factors once - eliminates colnumber² redundant conversions
  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)
  repli <- nlevels(replication)
  genotype <- nlevels(genotypes)
  
  # Pre-compute constants - avoids recalculating in nested loops
  CF_denom <- repli * genotype
  DFR <- repli - 1
  DFG <- genotype - 1
  DFE <- DFR * DFG
  repli_inv <- 1 / repli
  genotype_inv <- 1 / genotype
  
  # Convert to integer indices - rowsum() is faster with integer grouping
  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  
  # Pre-allocate result matrix - avoids growing vector (memory reallocation)
  phenotypic.cov <- matrix(0, nrow = colnumber, ncol = colnumber,
                           dimnames = list(headings, headings))
  
  # Vectorized computation using matrix operations
  for (i in seq_len(colnumber)) {
    trait1 <- data_mat[, i]
    
    # Pre-compute trait1 summaries for reuse in inner loop
    # rowsum() is 5-10x faster than tapply() - optimized C implementation
    sumch1 <- rowsum(trait1, gen_idx, reorder = FALSE)
    sumr1 <- rowsum(trait1, rep_idx, reorder = FALSE)
    GT1 <- sum(trait1)
    
    for (j in seq_len(colnumber)) {
      trait2 <- data_mat[, j]
      
      # rowsum() avoids tapply's S3 dispatch and list overhead
      sumch2 <- rowsum(trait2, gen_idx, reorder = FALSE)
      sumr2 <- rowsum(trait2, rep_idx, reorder = FALSE)
      GT2 <- sum(trait2)
      
      # Vectorized calculations
      CF <- (GT1 * GT2) / CF_denom
      
      # crossprod() is faster than sum(x * y) - direct BLAS call, no intermediate vector
      TSP <- crossprod(trait1, trait2)[1] - CF
      GSP <- crossprod(sumch1, sumch2)[1] * repli_inv - CF
      RSP <- crossprod(sumr1, sumr2)[1] * genotype_inv - CF
      
      ESP <- TSP - GSP - RSP
      
      # Compute mean products
      EMP <- ESP / DFE
      GMP <- GSP / DFG
      
      # Genetic and phenotypic covariance
      # No intermediate rounding preserves numerical precision
      GCov <- (GMP - EMP) * repli_inv
      phenotypic.cov[i, j] <- GCov + EMP
    }
  }
  
  return(phenotypic.cov)
}
