#' Convert dataframe to matrix
#'
#' @param data dataframe of weight
#'
#' @return A matrix
#' @export
#'
#' @examples
#' weight.mat(data = weight)
weight.mat<- function(data){
  datam<- data[,-1]
  m<- as.matrix(datam)
  return(m)
}
