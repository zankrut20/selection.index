#' Analysis of selection score and Genotype/Progeny ranking based on selection score
#'
#' @param data traits to be analyzed
#' @param bmat Selection/discriminant coefficient matrix
#' @param genotype Vector containing Genotypes/Treatments
#'
#' @return A data frame containing the genotypes with selection score
#' @importFrom stats aggregate
#' @export
#'
#' @examples
#' sel.score.rank(data = seldata[,3], bmat = 0.6455, genotype = seldata[,2])
sel.score.rank<- function(data, bmat, genotype){
  b<- matrix(bmat, ncol = 1)
  mean.df<- function(data, genotype){
    odr<- unique(genotype)
    geno<- list(factor(genotype, ordered(odr)))
    mean<- aggregate(data, geno, mean)
    matrix<- as.matrix(mean[,-1])
    return(matrix)
  }
  df<- mean.df(data, genotype)
  odr<- unique(genotype)
  sel.score<- df %*% b
  rank<- as.numeric(rank(as.vector(-sel.score)))
  rank.df<- data.frame("Genotype" = odr,
                       "Selection score" = sel.score,
                       "Rank" = rank)
  return(rank.df)
}
