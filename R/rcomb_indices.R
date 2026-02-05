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
  
  # Call C++ iterator to process all combinations
  result <- cpp_comb_iterator(
    pmat = pmat,
    gmat = gmat,
    wmat = wmat,
    comb_matrix = comb,
    wcol = wcol - 1,  # Convert to 0-indexed for C++
    const_factor = const_factor,
    PRE_constant = PRE_constant
  )
  
  # Construct data frame from C++ results
  df <- data.frame(
    ID = result$IDs,
    result$b_matrix,
    GA = result$GAs,
    PRE = result$PREs,
    Rank = rank(-result$PREs, ties.method = "min"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  return(df)
}
