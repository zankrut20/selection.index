#' Ranking of selection index
#'
#' @param list A list containing the selection indices components
#' @param i Number of top index
#'
#' @return A DataFrame containing the i ranked indecies
#' @export
#'
#' @examples
rank.index<- function(list, i){
  i = as.numeric(i)
  df<- do.call(rbind.data.frame, list)
  df$Rank<- as.numeric(rank(as.vector(-df$PRE)))
  df<- df[order(df$PRE, decreasing = TRUE),]
  return(df[1:i,])
}
