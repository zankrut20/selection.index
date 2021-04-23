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
#' gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' wmat<- weight.mat(weight)
#' comb.indices(ncomb = 1, pmat = pmat, gmat = gmat, wmat = wmat, w_col = 1, GAY = 1.075)
#'
comb.indices<- function(ncomb, pmat, gmat, wmat, wcol = 1, GAY){
  selection.index<- function(ID, pmat, gmat, wmat, GA){
    ID = toString(ID)
    p<- as.matrix(pmat)
    g<- as.matrix(gmat)
    w<- as.matrix(wmat)
    bmat<- solve(pmat) %*% gmat %*% wmat
    G<- 2.063 * t(bmat) %*% gmat %*% wmat / (t(bmat) %*% pmat %*% bmat)^0.5
    PRE<- if(missing(GA)){
      (G/G) * 100
    } else {
      (G/GA) * 100
    }
    result<- list("ID" = ID, "b" = matrix(round(bmat,4), nrow = 1),
                  "GA" = round(G,4), "PRE" = round(PRE,4))
    return(data.frame(result))
  }
  ncolmn<- ncol(pmat)
  comb<- t(combn(ncolmn, ncomb))
  indices<- list()
  for (i in 1:nrow(comb)) {
    as.numeric(ID<- paste0(comb[i,]))
    indices[[i]]<-selection.index(ID,
                                  pmat = pmat[comb[i,], comb[i,]],
                                  gmat = gmat[comb[i,], comb[i,]],
                                  wmat = wmat[comb[i,], wcol], GA = GAY)
  }
  df<- do.call(rbind.data.frame, indices)
  df$Rank<- as.numeric(rank(as.vector(-df$PRE)))
  # df<- df[order(df$PRE, decreasing = TRUE),]
  return(df)
}
