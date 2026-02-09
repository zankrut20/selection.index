#' Remove trait or trait combination from possible trait combinations of possible
#' Trait combinations
#'
#' @param ncomb Number of character combination
#' @param i remove trait or trait combination
#' @param pmat Phenotypic Variance Covariance Matrix
#' @param gmat Genotypic Variance Covariance Matrix
#' @param wmat Weight Matrix
#' @param wcol Respective weight column number of Weight Matrix
#' @param GAY Genetic Advance/Genetic Gain of base selection index
#'
#' @return Data frame of possible selection indices with per cent relative efficiency and ranking
#' @export
#'
#' @examples
#' gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' rcomb_indices(ncomb = 2, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,2:3], wcol = 1)
rcomb_indices<- function (ncomb, i, pmat, gmat, wmat, wcol = 1, GAY)
{
  # Convert matrices once
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  wmat <- as.matrix(wmat)
  
  ncolmn <- ncol(pmat)
  
  # Generate and filter combinations
  comb_all <- combn(ncolmn, ncomb)
  
  # Vectorized filtering using colSums
  keep_mask <- colSums(comb_all != i) == ncomb
  comb <- comb_all[, keep_mask, drop = FALSE]
  
  # Handle edge case: no combinations after filtering
  if (ncol(comb) == 0) {
    return(data.frame(ID = character(0), GA = numeric(0), 
                     PRE = numeric(0), Rank = numeric(0)))
  }
  
  # Pre-compute constants
  const_factor <- 2.063
  PRE_constant <- if(missing(GAY)) 100 else 100 / GAY
  ncomb_total <- ncol(comb)
  
  # Pre-allocate result storage
  IDs <- character(ncomb_total)
  b_list <- vector("list", ncomb_total)
  GAs <- numeric(ncomb_total)
  PREs <- numeric(ncomb_total)
  
  # Process each combination using math primitives
  for (j in seq_len(ncomb_total)) {
    # Get trait indices for this combination (1-indexed)
    idx <- comb[, j]
    IDs[j] <- paste(idx, collapse = ", ")
    
    # Extract submatrices for this combination
    P_sub <- cpp_extract_submatrix(pmat, idx)
    G_sub <- cpp_extract_submatrix(gmat, idx)
    w_sub <- cpp_extract_vector(wmat, idx, wcol - 1L)  # Convert to 0-indexed
    
    # Calculate selection index coefficients: b = P^(-1) * G * w
    Gw <- G_sub %*% w_sub
    b <- cpp_symmetric_solve(P_sub, Gw)
    
    # Calculate genetic advance: GA = const_factor * (b' * G * w) / sqrt(b' * P * b)
    numerator <- const_factor * cpp_quadratic_form(b, G_sub, w_sub)
    denominator <- sqrt(cpp_quadratic_form_sym(b, P_sub))
    GA <- numerator / denominator
    
    # Calculate percent relative efficiency
    PRE <- GA * PRE_constant
    
    # Store results (rounded to 4 decimals)
    b_list[[j]] <- round(as.vector(b), 4)
    GAs[j] <- round(GA, 4)
    PREs[j] <- round(PRE, 4)
  }
  
  # Convert b_list to matrix (pad with NA for shorter vectors)
  max_b_cols <- max(sapply(b_list, length))
  b_matrix <- matrix(NA_real_, nrow = ncomb_total, ncol = max_b_cols)
  colnames(b_matrix) <- paste0("b.", seq_len(max_b_cols))
  
  for (j in seq_len(ncomb_total)) {
    b_len <- length(b_list[[j]])
    b_matrix[j, 1:b_len] <- b_list[[j]]
  }
  
  # Construct result data frame
  df <- data.frame(
    ID = IDs,
    b_matrix,
    GA = GAs,
    PRE = PREs,
    Rank = rank(-PREs, ties.method = "min"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  return(df)
}
