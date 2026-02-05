#' @title Construction of selection indices based on number of character grouping
#'
#' @param ncomb Number of Characters/Traits group
#' @param pmat Phenotypic Variance-Covariance Matrix
#' @param gmat Genotypic Variance-Covariance Matrix
#' @param wmat Weight Matrix
#' @param wcol Weight column number incase more than one weights, by default its 1
#' @param GAY Genetic Advance of comparative Character/Trait i.e. Yield (Optional argument)
#' @return Data frame of all possible selection indices
#' @export
#' @importFrom utils combn
#' @examples
#' gmat<- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat<- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' wmat<- weight_mat(weight)
#' comb_indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1, GAY = 1.075)
#'
comb_indices<- function(ncomb, pmat, gmat, wmat, wcol = 1, GAY){
  # Convert matrices once
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  wmat <- as.matrix(wmat)
  
  # Generate combinations (keep as columns)
  ncolmn <- ncol(pmat)
  comb <- combn(ncolmn, ncomb)
  
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
