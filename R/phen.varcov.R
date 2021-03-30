#' Phenotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments
#' @param replication vector containing replication
#'
#' @return A Phenotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' phen.varcov(data=seldata[,3:9], genotypes=seldata$treat,replication=seldata$rep)
phen.varcov<- function (data, genotypes, replication)
{
  convrt<- function(data1) {
    data1<- as.data.frame(sapply(data1, as.numeric))
    data1<- as.list(data1)
    return(data1)
  }
  datam<- convrt(data)
  colnumber<- ncol(data)
  headings<- names(data)
  analysis<- function(genotypes, replication, trait1, trait2) {
    genotypes<- as.factor(genotypes)
    replication<- as.factor(replication)
    sumch1<- tapply(trait1, genotypes, sum)
    sumch2<- tapply(trait2, genotypes, sum)
    sumr1<- tapply(trait1, replication, sum)
    sumr2<- tapply(trait2, replication, sum)
    repli<- nlevels(replication)
    genotype<- nlevels(genotypes)
    GT1<- sum(trait1)
    GT2<- sum(trait2)
    CF<- (GT1 * GT2)/(repli * genotype)
    TSP<- round(sum(trait1 * trait2) - CF, 4)
    GSP<- round((sum(sumch1 * sumch2)/repli) - CF, 4)
    RSP<- round((sum(sumr1 * sumr2)/genotype) - CF, 4)
    ESP<- TSP - GSP - RSP
    DFR<- repli - 1
    DFG<- genotype - 1
    DFE<- DFR * DFG
    RMP<- round(RSP/DFR, 4)
    GMP<- round(GSP/DFG, 4)
    EMP<- round(ESP/DFE, 4)
    ECov<- round((GMP - EMP)/repli, 4)
    GCov<- round((GMP - EMP)/repli, 4)
    PCov<- GCov + ECov
    return(PCov)
  }
  phenotypic.cov<- c()
  index = 0
  for (i in 1:(colnumber)) {
    for (j in 1:colnumber) {
      index = index + 1
      phenotypic.cov[index]<- analysis(genotypes, replication,
                                       datam[[i]], datam[[j]])
    }
  }
  matrix1<- matrix(phenotypic.cov, nrow = colnumber,
                   dimnames = list(headings, headings))
  return(matrix1)
}
