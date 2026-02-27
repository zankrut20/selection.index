#' Convert dataframe to matrix
#'
#' @param data dataframe of weight
#'
#' @return A matrix
#' @export
#'
#' @examples
#' weight_mat(data = weight)
weight_mat <- function(data) {
  m <- as.matrix(data[, -1, drop = FALSE])
  m
}
