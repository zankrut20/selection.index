#' Genetic Advance for PRE
#'
#' @param phen_mat phenotypic matrix value of desired characters
#' @param gen_mat genotypic matrix value of desired characters
#' @param weight_mat weight matrix value of desired characters
#' @return Genetic advance of character or character combinations
#' @export
#'
#' @examples
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' gen_advance(phen_mat = pmat[1, 1], gen_mat = gmat[1, 1], weight_mat = weight[1, 2])
gen_advance <- function(phen_mat, gen_mat, weight_mat) {
  p <- as.matrix(phen_mat)
  g <- as.matrix(gen_mat)
  w <- as.matrix(weight_mat)

  bmat <- solve(p, g %*% w)

  numerator <- 2.063 * crossprod(bmat, g %*% w)
  denominator <- sqrt(crossprod(bmat, p %*% bmat))

  GA <- round(numerator / denominator, 4)

  GA
}
