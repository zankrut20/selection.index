#' Convert dataframe to matrix
#'
#' @param data dataframe of weight
#'
#' @return A matrix
#' @export
#'
#' @examples
weight.mat<- function(data){
  rownames(data)<- paste(data[,1])
  datam<- data[,-1]
  m<- as.matrix(datam)
  return(m)
}
