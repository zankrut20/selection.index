#' Experimental Design Statistics Engine
#'
#' @description
#' Modular engine for experimental design statistics supporting RCBD and Latin Square designs.
#' Computes correction factors, sums of products, mean products, and degrees of freedom
#' for variance-covariance analysis and ANOVA statistics.
#'
#' This function eliminates code duplication across gen.varcov(), phen.varcov(), and 
#' mean.performance() by providing a centralized, optimized implementation.
#'
#' @param trait1 Numeric vector of first trait observations
#' @param trait2 Numeric vector of second trait observations (default: trait1 for variance)
#' @param genotypes Integer vector of genotype/treatment indices
#' @param replications Integer vector of replication/block indices (for RCBD) or row indices (for LSD)
#' @param columns Integer vector of column indices (required for Latin Square Design only)
#' @param design_type Character string specifying design type: "RCBD" (default) or "LSD" (Latin Square)
#' @param calc_type Character string specifying calculation type:
#'   \itemize{
#'     \item \code{"sums_of_products"} - Returns CF, TSP, GSP, RSP, ESP (for covariance)
#'     \item \code{"mean_products"} - Returns GMP, EMP (for variance components)
#'     \item \code{"all"} - Returns all components (default)
#'     \item \code{"anova_stats"} - Returns degrees of freedom and mean squares
#'   }
#'
#' @return List containing RCBD statistics based on calc_type:
#'   \itemize{
#'     \item \code{CF} - Correction factor: (GT1 × GT2) / (r × t)
#'     \item \code{TSP} - Total sum of products
#'     \item \code{GSP} - Genotype sum of products  
#'     \item \code{RSP} - Replication sum of products
#'     \item \code{ESP} - Error sum of products
#'     \item \code{GMP} - Genotype mean product
#'     \item \code{EMP} - Error mean product
#'     \item \code{DFG} - Degrees of freedom for genotypes (t - 1)
#'     \item \code{DFR} - Degrees of freedom for replications (r - 1)
#'     \item \code{DFE} - Degrees of freedom for error (r - 1)(t - 1)
#'     \item \code{n_genotypes} - Number of genotypes
#'     \item \code{n_replications} - Number of replications
#'   }
#'
#' @details
#' The function uses optimized Base R operations:
#' \itemize{
#'   \item \code{rowsum()} for fast grouped summations (5-10x faster than tapply)
#'   \item \code{crossprod()} for efficient matrix products (faster than sum(x*y))
#'   \item Pre-computed constants to avoid repeated calculations
#' }
#'
#' **RCBD Model:** Y_ij = μ + τ_i + β_j + ε_ij  
#' where τ_i = genotype effect, β_j = block effect, ε_ij = error
#' 
#' **LSD Model:** Y_ijk = μ + τ_i + ρ_j + γ_k + ε_ijk  
#' where τ_i = genotype effect, ρ_j = row effect, γ_k = column effect, ε_ijk = error
#'
#' @references
#' Cochran, W. G., & Cox, G. M. (1957). Experimental designs (2nd ed.). Wiley.
#' 
#' Steel, R. G. D., & Torrie, J. H. (1980). Principles and procedures of statistics: 
#' A biometrical approach (2nd ed.). McGraw-Hill.
#'
#' @keywords internal
design.stats <- function(trait1, trait2 = trait1, genotypes, replications, 
                        columns = NULL,
                        design_type = c("RCBD", "LSD"),
                        calc_type = c("all", "sums_of_products", "mean_products", "anova_stats")) {
  
  design_type <- match.arg(design_type)
  calc_type <- match.arg(calc_type)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  
  # OPTIMIZATION: Ensure numeric vectors (handle factors/characters)
  # storage.mode assignment is faster than as.numeric() for already-numeric data
  if (!is.numeric(trait1)) trait1 <- as.numeric(trait1)
  if (!is.numeric(trait2)) trait2 <- as.numeric(trait2)
  storage.mode(trait1) <- "numeric"
  storage.mode(trait2) <- "numeric"
  
  # OPTIMIZATION: Count levels once
  # Using unique() on integer vectors is faster than nlevels(factor)
  n_genotypes <- length(unique(genotypes))
  n_replications <- length(unique(replications))
  
  if (design_type == "RCBD") {
    # ========== RCBD CALCULATIONS ==========
    # OPTIMIZATION: Pre-compute constants
    # Avoid repeated arithmetic - division is ~10x slower than multiplication
    n_obs <- n_genotypes * n_replications
    repli_inv <- 1 / n_replications
    genotype_inv <- 1 / n_genotypes
    
    # Degrees of freedom
    DFG <- n_genotypes - 1
    DFR <- n_replications - 1
    DFE <- DFG * DFR
    
    # OPTIMIZATION: Use rowsum() for grouped sums (5-10x faster than tapply)
    # rowsum() is .Internal primitive optimized in C
    sumch1 <- rowsum(trait1, genotypes, reorder = FALSE)
    sumch2 <- rowsum(trait2, genotypes, reorder = FALSE)
    sumr1 <- rowsum(trait1, replications, reorder = FALSE)
    sumr2 <- rowsum(trait2, replications, reorder = FALSE)
    
    # Grand totals
    GT1 <- sum(trait1)
    GT2 <- sum(trait2)
    
    # Correction Factor
    CF <- (GT1 * GT2) / n_obs
    
    # Return early for anova_stats (no need to compute products)
    if (calc_type == "anova_stats") {
      return(list(
        DFG = DFG,
        DFR = DFR,
        DFE = DFE,
        n_genotypes = n_genotypes,
        n_replications = n_replications,
        CF = CF,
        design_type = "RCBD"
      ))
    }
    
    # OPTIMIZATION: Use crossprod() for sum of products
    # Faster than sum(x * y) - direct BLAS call, no intermediate vector
    TSP <- crossprod(trait1, trait2)[1] - CF
    GSP <- crossprod(sumch1, sumch2)[1] * repli_inv - CF
    RSP <- crossprod(sumr1, sumr2)[1] * genotype_inv - CF
    
    # Error sum of products
    ESP <- TSP - GSP - RSP
    
    # Return early for sums_of_products
    if (calc_type == "sums_of_products") {
      return(list(
        CF = CF,
        TSP = TSP,
        GSP = GSP,
        RSP = RSP,
        ESP = ESP,
        DFG = DFG,
        DFR = DFR,
        DFE = DFE,
        n_genotypes = n_genotypes,
        n_replications = n_replications,
        design_type = "RCBD"
      ))
    }
    
    # Mean products
    GMP <- GSP / DFG
    EMP <- ESP / DFE
    
    # Return for mean_products
    if (calc_type == "mean_products") {
      return(list(
        GMP = GMP,
        EMP = EMP,
        DFG = DFG,
        DFR = DFR,
        DFE = DFE,
        n_genotypes = n_genotypes,
        n_replications = n_replications,
        design_type = "RCBD"
      ))
    }
    
    # Return all (default)
    return(list(
      CF = CF,
      TSP = TSP,
      GSP = GSP,
      RSP = RSP,
      ESP = ESP,
      GMP = GMP,
      EMP = EMP,
      DFG = DFG,
      DFR = DFR,
      DFE = DFE,
      n_genotypes = n_genotypes,
      n_replications = n_replications,
      design_type = "RCBD"
    ))
    
  } else {
    # ========== LATIN SQUARE DESIGN CALCULATIONS ==========
    n_columns <- length(unique(columns))
    
    # For LSD: n_replications represents rows, n_columns represents columns
    # t = number of treatments (genotypes)
    t <- n_genotypes
    
    # OPTIMIZATION: Pre-compute constants
    n_obs <- t * t  # In LSD, typically t×t observations
    t_inv <- 1 / t
    
    # Degrees of freedom for LSD
    DFG <- t - 1           # Treatments (genotypes)
    DFR <- t - 1           # Rows
    DFC <- t - 1           # Columns
    DFE <- (t - 1) * (t - 2)  # Error: (t-1)(t-2)
    
    # OPTIMIZATION: Use rowsum() for grouped sums
    sumch1 <- rowsum(trait1, genotypes, reorder = FALSE)
    sumch2 <- rowsum(trait2, genotypes, reorder = FALSE)
    sumr1 <- rowsum(trait1, replications, reorder = FALSE)  # rows
    sumr2 <- rowsum(trait2, replications, reorder = FALSE)
    sumc1 <- rowsum(trait1, columns, reorder = FALSE)
    sumc2 <- rowsum(trait2, columns, reorder = FALSE)
    
    # Grand totals
    GT1 <- sum(trait1)
    GT2 <- sum(trait2)
    
    # Correction Factor
    CF <- (GT1 * GT2) / n_obs
    
    # Return early for anova_stats
    if (calc_type == "anova_stats") {
      return(list(
        DFG = DFG,
        DFR = DFR,
        DFC = DFC,
        DFE = DFE,
        n_genotypes = t,
        n_rows = t,
        n_columns = t,
        CF = CF,
        design_type = "LSD"
      ))
    }
    
    # OPTIMIZATION: Use crossprod() for sum of products
    TSP <- crossprod(trait1, trait2)[1] - CF
    GSP <- crossprod(sumch1, sumch2)[1] * t_inv - CF  # Genotype/Treatment SP
    RSP <- crossprod(sumr1, sumr2)[1] * t_inv - CF    # Row SP
    CSP <- crossprod(sumc1, sumc2)[1] * t_inv - CF    # Column SP
    
    # Error sum of products for LSD
    ESP <- TSP - GSP - RSP - CSP
    
    # Return early for sums_of_products
    if (calc_type == "sums_of_products") {
      return(list(
        CF = CF,
        TSP = TSP,
        GSP = GSP,
        RSP = RSP,
        CSP = CSP,
        ESP = ESP,
        DFG = DFG,
        DFR = DFR,
        DFC = DFC,
        DFE = DFE,
        n_genotypes = t,
        n_rows = t,
        n_columns = t,
        design_type = "LSD"
      ))
    }
    
    # Mean products
    GMP <- GSP / DFG
    EMP <- ESP / DFE
    
    # Return for mean_products
    if (calc_type == "mean_products") {
      return(list(
        GMP = GMP,
        EMP = EMP,
        DFG = DFG,
        DFR = DFR,
        DFC = DFC,
        DFE = DFE,
        n_genotypes = t,
        n_rows = t,
        n_columns = t,
        design_type = "LSD"
      ))
    }
    
    # Return all (default)
    return(list(
      CF = CF,
      TSP = TSP,
      GSP = GSP,
      RSP = RSP,
      CSP = CSP,
      ESP = ESP,
      GMP = GMP,
      EMP = EMP,
      DFG = DFG,
      DFR = DFR,
      DFC = DFC,
      DFE = DFE,
      n_genotypes = t,
      n_rows = t,
      n_columns = t,
      design_type = "LSD"
    ))
  }
}
