#' Construction of Selection Index
#'
#' @param ID Name of Selection Index
#' @param phen_mat Phenotypic Variance-Covariance Matrix
#' @param gen_mat Genotypic Variance-Covariance Matrix
#' @param weight_mat Weight Matrix
#' @param GAY Genetic Advance of Yield
#'
#' @return List of Selection Index Components
#' @export
#'
#' @examples
selection.index<- function(ID, phen_mat, gen_mat, weight_mat, GAY)
{
  ID = toString(ID)
  p<- as.matrix(phen_mat)
  g<- as.matrix(gen_mat)
  w<- as.matrix(weight_mat)
  bmat<- solve(phen_mat) %*% gen_mat %*% weight_mat
  GA<- 2.063 * t(bmat) %*% g %*% w / (t(bmat) %*% p %*% bmat)^0.5
  PRE<- if(missing(GAY)){
    (GA/GA) * 100
  } else {
    (GA/GAY) * 100
  }
  result<- list("ID" = ID, "b" = matrix(round(bmat,4), nrow = 1),
                "GA" = round(GA,4), "PRE" = round(PRE,4))
  return(result)
}
