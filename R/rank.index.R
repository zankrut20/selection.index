#' Ranking of selection index
#'
#' @param list A list containing the selection indices components
#' @param i Number of top index
#'
#' @return A DataFrame containing the i ranked indecies
#' @export
#'
#' @examples
#' gmat<- gen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat<- phen.varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' si<- list()
#' si[[1]]<- sel.index(ID = 1, phen_mat = pmat[1,1], gen_mat = gmat[1,1], weight_mat = weight[1,2])
#' si[[2]]<- sel.index(ID = 1, phen_mat = pmat[2,2], gen_mat = gmat[2,2], weight_mat = weight[2,2])
#' si[[3]]<- sel.index(ID = 1, phen_mat = pmat[3,3], gen_mat = gmat[3,3], weight_mat = weight[3,2])
#' si[[4]]<- sel.index(ID = 1, phen_mat = pmat[4,4], gen_mat = gmat[4,4], weight_mat = weight[4,2])
#' rank.index(list = si, i = 2)
rank.index<- function(list, i){
  i = as.numeric(i)
  df<- do.call(rbind.data.frame, list)
  df$Rank<- as.numeric(rank(as.vector(-df$PRE)))
  df<- df[order(df$PRE, decreasing = TRUE),]
  return(df[1:i,])
}
