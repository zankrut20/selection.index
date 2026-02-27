#' Synthetic Maize Phenotypic Data
#'
#' A synthetic dataset containing multi-environment phenotypic records for 100
#' maize genotypes. Designed to demonstrate phenotypic selection index
#' calculations and variance-covariance matrix estimations.
#'
#' @name maize_pheno
#' @docType data
#' @format A data frame with 600 rows and 6 variables:
#' \describe{
#'   \item{Genotype}{Factor representing 100 unique maize genotypes.}
#'   \item{Environment}{Factor representing 2 distinct growing environments.}
#'   \item{Block}{Factor representing 3 replicate blocks within each environment.}
#'   \item{Yield}{Numeric vector of grain yield in kg/ha.}
#'   \item{PlantHeight}{Numeric vector of plant height in centimeters.}
#'   \item{DaysToMaturity}{Numeric vector of days to physiological maturity.}
#' }
#' @source Simulated for the selection.index package to provide reproducible examples.
#' @usage data(maize_pheno)
#' @keywords datasets
#' @examples
#' data(maize_pheno)
#' head(maize_pheno)
NULL

#' Synthetic Maize Genomic Data
#'
#' A synthetic dataset containing Single Nucleotide Polymorphism (SNP) marker
#' data for the 100 maize genotypes found in maize_pheno. Designed for
#' testing genomic selection indices (GESIM/RGESIM) and relationship matrices.
#'
#' @name maize_geno
#' @docType data
#' @format A numeric matrix with 100 rows (genotypes) and 500 columns (SNP markers):
#' \describe{
#'   \item{Rows}{100 genotypes, matching the Genotype column in maize_pheno.}
#'   \item{Columns}{500 simulated SNP markers, coded as 0, 1, or 2.}
#' }
#' @source Simulated for the selection.index package to provide reproducible examples.
#' @usage data(maize_geno)
#' @keywords datasets
#' @examples
#' data(maize_geno)
#' dim(maize_geno)
NULL

#' Selection Index DataSet
#'
#' A dataset containing the data of three replications and 48 progenies
#' with 7 different traits.
#'
#' \itemize{
#' \item rep. Replications
#' \item treat. Treatments/Genotypes
#' \item sypp. Seed Yield per Plant
#' \item dtf. Days to 50% Flowering
#' \item rpp. Racemes per Plant
#' \item ppr. Pods per Raceme
#' \item ppp. Pods per Plant
#' \item spp. Seeds per Pod
#' \item pw. Pods Weight
#' }
#' @docType data
#' @name seldata
#' @usage data(seldata)
#' @format A data frame with 75 rows and 9 columns
NULL

#' Weight dataset
#'
#' A dataset containing the data of 2 different weights namely euqal weight and
#' broad sense heritability
#'
#' \itemize{
#' \item EW. Equal Weight
#' \item h2. Broad Sense Heritability
#' }
#' @docType data
#' @name weight
#' @usage data(weight)
#' @format A data frame with 7 rows and 3 columns
NULL
