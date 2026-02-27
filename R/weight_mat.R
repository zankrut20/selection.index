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
  # OPTIMIZATION: Single-step subsetting and conversion
  # Avoids: Creating intermediate variable 'datam' (extra memory allocation)
  # Why faster: as.matrix(data[,-1]) does subsetting and conversion in one pass
  # Note: Column subsetting with negative index is already efficient in base R
  m <- as.matrix(data[, -1, drop = FALSE])
  m
}
