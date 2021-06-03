#' Remove trait or trait combination from possible trait combinations of possible
#' Trait combinations
#'
#' @param ncomb Number of character combination
#' @param i remove trait or trait combination
#' @param pmat Phenotypic Variance Covariance Matrix
#' @param gmat Genotypic Variane Covariance Matrix
#' @param wmat Weight Matrix
#' @param wcol Respective weight column number of Weight Matrix
#' @param GAY Genetic Advance/Genetic Gain of base selection index
#'
#' @return Data frame of possible selection indices with per cent relative efficiency and ranking
#' @export
#'
#' @examples
#' gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' rcomb.indices(ncomb = 2, i = 1, pmat = pmat, gmat = gmat, wmat = weight[,2:3], wcol = 1)
rcomb.indices<- function (ncomb, i, pmat, gmat, wmat, wcol = 1, GAY)
{
  selection.index <- function(ID, pmat, gmat, wmat, GA) {
    ID = toString(ID)
    p <- as.matrix(pmat)
    g <- as.matrix(gmat)
    w <- as.matrix(wmat)
    bmat <- solve(pmat) %*% gmat %*% wmat
    G <- 2.063 * t(bmat) %*% gmat %*% wmat/(t(bmat) %*%
                                              pmat %*% bmat)^0.5
    PRE <- if (missing(GA)) {
      (G/G) * 100
    }
    else {
      (G/GA) * 100
    }
    result <- list(ID = ID, b = matrix(round(bmat, 4), nrow = 1),
                   GA = round(G, 4), PRE = round(PRE, 4))
    return(data.frame(result))
  }
  ncolmn <- ncol(pmat)
  trait.remove<- function(n, c, i){
    a<- t(combn(n,c))
    b<- if (c == 1){
      matrix(a[rowSums(a !=i) == c, ])
    } else if ( n - c > 1) {
      a[rowSums(a !=i) == c, ]
    } else {
      t(matrix(a[rowSums(a !=i) == c, ]))
    }
    return(b)
  }
  comb <- trait.remove(ncolmn, ncomb, i)
  indices <- list()
  for (i in 1:nrow(comb)) {
    as.numeric(ID <- paste0(comb[i, ]))
    indices[[i]] <- selection.index(ID, pmat = pmat[comb[i,
    ], comb[i, ]], gmat = gmat[comb[i, ], comb[i, ]],
    wmat = wmat[comb[i, ], wcol], GA = GAY)
  }
  df <- do.call(rbind.data.frame, indices)
  return(df)
}
