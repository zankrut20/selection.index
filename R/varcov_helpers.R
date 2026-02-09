#' Calculate Variance-Covariance Components (Internal Helper)
#' 
#' @description 
#' Internal function to compute genotypic or phenotypic variance-covariance
#' matrices using math primitives. Used by gen_varcov() and phen_varcov().
#' 
#' @param data_mat Numeric matrix of trait data
#' @param gen_idx Integer vector of genotype indices
#' @param rep_idx Integer vector of replication indices
#' @param col_idx Integer vector of column indices (for LSD, optional)
#' @param main_idx Integer vector of main plot indices (for SPD, optional)
#' @param design_type Integer design code: 1=RCBD, 2=LSD, 3=SPD
#' @param cov_type Integer: 1=genotypic, 2=phenotypic
#' 
#' @return Symmetric variance-covariance matrix
#' 
#' @keywords internal
#' @noRd
.calculate_varcov <- function(data_mat, gen_idx, rep_idx,
                               col_idx = NULL, main_idx = NULL,
                               design_type = 1L, cov_type = 1L) {
  
  n_traits <- ncol(data_mat)
  n_obs <- nrow(data_mat)
  
  if (cov_type == 2L) {
    # Phenotypic variance-covariance: simple covariance matrix
    # But need to adjust for experimental design
    
    # Get ANOVA components
    anova_result <- .calculate_anova(data_mat, gen_idx, rep_idx,
                                     col_idx, main_idx, design_type)
    
    MSG <- anova_result$MSG
    MSE <- anova_result$MSE
    
    # Phenotypic variance = MSG (contains both genetic and error variance)
    # This is the among-genotype variance component
    return(MSG)
    
  } else {
    # Genotypic variance-covariance
    
    # Get ANOVA components
    anova_result <- .calculate_anova(data_mat, gen_idx, rep_idx,
                                     col_idx, main_idx, design_type)
    
    MSG <- anova_result$MSG
    MSE <- anova_result$MSE
    
    # Calculate genotypic variance-covariance based on design
    if (design_type == 1L) {  # RCBD
      n_rep <- length(unique(rep_idx))
      Vg <- (MSG - MSE) / n_rep
      
    } else if (design_type == 2L) {  # LSD
      n_rep <- length(unique(rep_idx))
      Vg <- (MSG - MSE) / n_rep
      
    } else if (design_type == 3L) {  # SPD
      n_rep <- length(unique(rep_idx))
      n_main <- length(unique(main_idx))
      Vg <- (MSG - MSE) / (n_rep * n_main)
    }
    
    return(Vg)
  }
}
