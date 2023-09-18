#' @title Mean performance of phenotypic data
#'
#' @param data data for analysis
#' @param genotypes genotypes vector
#' @param replications replication vector
#'
#' @return Dataframe of mean performance analysis
#' @export
#' @importFrom stats aggregate anova lm qt
#' @examples
#' meanPerformance(data = seldata[, 3:9], genotypes = seldata[, 2], replications = seldata[, 1])
#'
#'

meanPerformance <- function(data, genotypes, replications){
  # Calculate average results from replicaions
  meanData <- function(data, genotypes){
    odr <- unique(genotypes)
    genotype <- list(factor(genotypes, ordered(odr)))
    mean <- aggregate(data, genotype, mean)
    matrix <- as.matrix(mean[,-1])
    matrix <- round(matrix, 4)
    result <- data.frame("Genotypes" = odr, matrix)
    return(data.frame(result))
  }
  # Calculate the performance
  performance <- function (data, treat, rep)
  {
    convert <- function(data) {
      data <- as.data.frame(sapply(data, as.numeric))
      data <- as.list(data)
      return(data)
    }
    analysis <- function(data, treat, rep) {
      data <- as.numeric(data)
      genotype <- as.factor(treat)
      replication <- as.factor(rep)
      r <- nlevels(replication)
      model <- lm(data ~ replication + genotype)
      anova.model <- anova(model)
      MinMax <- meanData(data = data, genotypes = treat)
      Maxi <- round(max(MinMax[, -1]), 4)
      Mini <- round(min(MinMax[, -1]), 4)
      rang <- paste(round(Maxi, 4), " - ", round(Mini, 4))
      GM <- round(mean(data), 4)
      EMS <- round(anova.model[3, 3],4)
      SD <- round(sqrt(EMS),4)
      SEm <- round(sqrt(EMS/r), 4)
      CV<- round(SD/GM*100, 4)
      CD5 <- round(sqrt(EMS/r) * sqrt(2) * abs(qt(0.025, anova.model[3, 1])), 4)
      if (anova.model[2, 5] > 0.05) {
        CD5 <- paste(round(sqrt(EMS/r) * sqrt(2) * abs(qt(0.025,
                                                          anova.model[3, 1])), 4), "NS")
      }
      CD1 <- round(sqrt(EMS/r) * sqrt(2) * abs(qt(0.005, anova.model[3,
                                                                     1])), 4)
      if (anova.model[2, 5] > 0.01) {
        CD1 <- paste(round(sqrt(EMS/r) * sqrt(2) * abs(qt(0.005,
                                                          anova.model[3, 1])), 4), "NS")
      }
      GV <- round((anova.model[2, 3] - EMS)/r, 4)
      PV <- round(GV + EMS, 4)
      hs <- round((GV/PV), 4)
      matri <- matrix(data = c(Mini, Maxi, GM, CV, SEm, CD5, CD1, hs, hs*100),
                      dimnames = list(c("Min", "Max", "Grand Mean", "CV %",
                                        "SEm",
                                        "CD 5%",
                                        "CD 1%",
                                        "Heritability",
                                        "Heritability in %")),
                      nrow = 9)
      table1 <- as.data.frame(matri, useNa = F)
      my.list <- list(table1)
      return(my.list)
    }
    fiftn <- convert(data)
    colnumber <- ncol(data)
    output <- list()
    for (j in 1:colnumber) {
      output[[j]] <- analysis(fiftn[[j]], treat,
                              rep)
    }
    names(output) <- names(data)
    return(output)
  }
  meandf <- meanData(data = data, genotypes = genotypes)
  performanceList <- performance(data = data, treat = genotypes, rep = replications)
  performanceBind <- cbind.data.frame(performanceList, row.names = NULL)
  Genotypes <- matrix(c("Min", "Max", "GM", "CV (%)", "SEm", "CD 5%", "CD 1%", "Heritability", "Heritability(%)"),
                      ncol = 1)
  FinalPerformanceBind <- cbind.data.frame(Genotypes, performanceBind, row.names = NULL)
  colnames(FinalPerformanceBind) <- colnames(meandf)
  result <- rbind.data.frame(meandf, FinalPerformanceBind)
  return(result)
}
