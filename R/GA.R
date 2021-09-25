#' Genetic Advance for PRE
#'
#' @param phen_mat phenotypic matrix value of desired characters
#' @param gen_mat genotypic matrix value of desired characters
#' @param weight_mat weight matrix value of desired characters
#'
#' @examples
#' gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' sel.index(phen_mat = pmat[1,1], gen_mat = gmat[1,1], weight_mat = weight[1,2])
gen.advance<- function(phen_mat, gen_mat, weight_mat)
{
  ID = toString(ID)
  p<- as.matrix(phen_mat)
  g<- as.matrix(gen_mat)
  w<- as.matrix(weight_mat)
  bmat<- solve(phen_mat) %*% gen_mat %*% weight_mat
  GA<- 2.063 * t(bmat) %*% g %*% w / (t(bmat) %*% p %*% bmat)^0.5
  return(as.numeric(GA))
}
