#' Experimental Design Statistics Engine
#'
#' @description
#' Modular engine for experimental design statistics supporting RCBD, Latin Square, and Split Plot designs.
#' Computes correction factors, sums of products, mean products, and degrees of freedom
#' for variance-covariance analysis and ANOVA statistics.
#'
#' This function eliminates code duplication across gen.varcov(), phen.varcov(), and 
#' mean.performance() by providing a centralized, optimized implementation.
#'
#' @param trait1 Numeric vector of first trait observations
#' @param trait2 Numeric vector of second trait observations (default: trait1 for variance)
#' @param genotypes Integer vector of genotype/treatment indices (sub-plot treatments in SPD)
#' @param replications Integer vector of replication/block indices (for RCBD) or row indices (for LSD)
#' @param columns Integer vector of column indices (required for Latin Square Design only)
#' @param main_plots Integer vector of main plot treatment indices (required for Split Plot Design only)
#' @param design_type Character string specifying design type: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot)
#' @param calc_type Character string specifying calculation type:
#'   \itemize{
#'     \item \code{"sums_of_products"} - Returns CF, TSP, GSP, RSP, ESP (for covariance)
#'     \item \code{"mean_products"} - Returns GMP, EMP (for variance components)
#'     \item \code{"all"} - Returns all components (default)
#'     \item \code{"anova_stats"} - Returns degrees of freedom and mean squares
#'   }
#'
#' @return List containing design statistics based on calc_type and design:
#'   \itemize{
#'     \item \code{CF} - Correction factor
#'     \item \code{TSP} - Total sum of products
#'     \item \code{GSP} - Genotype/Sub-plot sum of products  
#'     \item \code{RSP} - Replication sum of products
#'     \item \code{MSP} - Main plot sum of products (SPD only)
#'     \item \code{IMSP} - Main plot × Replication interaction SP (SPD only)
#'     \item \code{ESP} - Error sum of products (sub-plot error for SPD)
#'     \item \code{ESP_MAIN} - Main plot error sum of products (SPD only)
#'     \item \code{GMP} - Genotype mean product
#'     \item \code{EMP} - Error mean product (sub-plot error for SPD)
#'     \item \code{EMP_MAIN} - Main plot error mean product (SPD only)
#'     \item \code{DFG} - Degrees of freedom for genotypes/sub-plots
#'     \item \code{DFR} - Degrees of freedom for replications
#'     \item \code{DFM} - Degrees of freedom for main plots (SPD only)
#'     \item \code{DFE} - Degrees of freedom for error (sub-plot error for SPD)
#'     \item \code{DFE_MAIN} - Degrees of freedom for main plot error (SPD only)
#'     \item \code{n_genotypes} - Number of genotypes/sub-plot treatments
#'     \item \code{n_replications} - Number of replications
#'     \item \code{n_main_plots} - Number of main plot treatments (SPD only)
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
#' **SPD Model:** Y_ijk = μ + ρ_i + α_j + δ_ij + τ_k + (ατ)_jk + ε_ijk  
#' where ρ_i = block effect, α_j = main plot effect, δ_ij = main plot error,
#' τ_k = sub-plot effect (genotype), (ατ)_jk = interaction, ε_ijk = sub-plot error
#'
#' @references
#' Cochran, W. G., & Cox, G. M. (1957). Experimental designs (2nd ed.). Wiley.
#' 
#' Steel, R. G. D., & Torrie, J. H. (1980). Principles and procedures of statistics: 
#' A biometrical approach (2nd ed.). McGraw-Hill.
#' 
#' Gomez, K. A., & Gomez, A. A. (1984). Statistical procedures for agricultural research (2nd ed.). Wiley.
#'
#' @keywords internal
design_stats <- function(trait1, trait2 = trait1, genotypes, replications, 
                        columns = NULL, main_plots = NULL,
                        design_type = c("RCBD", "LSD", "SPD"),
                        calc_type = c("all", "sums_of_products", "mean_products", "anova_stats")) {
  
  design_type <- match.arg(design_type)
  calc_type <- match.arg(calc_type)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  
  # Validate Split Plot Design requirements
  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
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
    
  } else if (design_type == "LSD") {
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
  } else if (design_type == "SPD") {
    # ========== SPLIT PLOT DESIGN CALCULATIONS ==========
    # SPD Structure:
    # - Replications (blocks): r
    # - Main plot treatments: a
    # - Sub-plot treatments (genotypes): b
    # - Total observations: r × a × b
    
    # OPTIMIZATION: Count levels once
    n_main_plots <- length(unique(main_plots))
    n_sub_plots <- n_genotypes  # genotypes are sub-plot treatments
    n_obs <- length(trait1)
    
    # Pre-compute constants
    r <- n_replications
    a <- n_main_plots
    b <- n_sub_plots
    
    # Degrees of freedom for SPD
    DFR <- r - 1                      # Replications
    DFM <- a - 1                      # Main plot treatments
    DFE_MAIN <- DFR * DFM             # Main plot error (Block × Main plot interaction)
    DFG <- b - 1                      # Sub-plot treatments (genotypes)
    DFIM <- DFM * DFG                 # Main × Sub interaction
    DFE <- a * (b - 1) * DFR          # Sub-plot error
    
    # OPTIMIZATION: Use rowsum() for grouped sums
    # Grand totals
    GT1 <- sum(trait1)
    GT2 <- sum(trait2)
    
    # Correction Factor
    CF <- (GT1 * GT2) / n_obs
    
    # Sum by replications (blocks)
    sumr1 <- rowsum(trait1, replications, reorder = FALSE)
    sumr2 <- rowsum(trait2, replications, reorder = FALSE)
    
    # Sum by main plots
    summ1 <- rowsum(trait1, main_plots, reorder = FALSE)
    summ2 <- rowsum(trait2, main_plots, reorder = FALSE)
    
    # Sum by sub-plots (genotypes)
    sums1 <- rowsum(trait1, genotypes, reorder = FALSE)
    sums2 <- rowsum(trait2, genotypes, reorder = FALSE)
    
    # Sum by main plot within replications (for main plot error)
    # Create combined factor: replication × main_plot
    rep_main_factor <- paste(replications, main_plots, sep = "_")
    sumrm1 <- rowsum(trait1, rep_main_factor, reorder = FALSE)
    sumrm2 <- rowsum(trait2, rep_main_factor, reorder = FALSE)
    
    # Sum by genotype × main plot (for interaction)
    main_sub_factor <- paste(main_plots, genotypes, sep = "_")
    summs1 <- rowsum(trait1, main_sub_factor, reorder = FALSE)
    summs2 <- rowsum(trait2, main_sub_factor, reorder = FALSE)
    
    # Return early for anova_stats
    if (calc_type == "anova_stats") {
      return(list(
        DFR = DFR,
        DFM = DFM,
        DFE_MAIN = DFE_MAIN,
        DFG = DFG,
        DFIM = DFIM,
        DFE = DFE,
        n_replications = r,
        n_main_plots = a,
        n_genotypes = b,
        CF = CF,
        design_type = "SPD"
      ))
    }
    
    # OPTIMIZATION: Use crossprod() for sum of products
    TSP <- crossprod(trait1, trait2)[1] - CF
    
    # Replication sum of products
    RSP <- crossprod(sumr1, sumr2)[1] / (a * b) - CF
    
    # Main plot sum of products
    MSP <- crossprod(summ1, summ2)[1] / (r * b) - CF
    
    # Sub-plot (genotype) sum of products
    GSP <- crossprod(sums1, sums2)[1] / (r * a) - CF
    
    # Replication × Main plot sum of products (for main plot error)
    RMSP <- crossprod(sumrm1, sumrm2)[1] / b - CF
    
    # Main plot error sum of products
    ESP_MAIN <- RMSP - RSP - MSP
    
    # Main × Sub interaction sum of products
    IMSP <- crossprod(summs1, summs2)[1] / r - CF - MSP - GSP
    
    # Sub-plot error sum of products
    ESP <- TSP - RSP - MSP - ESP_MAIN - GSP - IMSP
    
    # Return early for sums_of_products
    if (calc_type == "sums_of_products") {
      return(list(
        CF = CF,
        TSP = TSP,
        RSP = RSP,
        MSP = MSP,
        GSP = GSP,
        IMSP = IMSP,
        ESP_MAIN = ESP_MAIN,
        ESP = ESP,
        DFR = DFR,
        DFM = DFM,
        DFE_MAIN = DFE_MAIN,
        DFG = DFG,
        DFIM = DFIM,
        DFE = DFE,
        n_replications = r,
        n_main_plots = a,
        n_genotypes = b,
        design_type = "SPD"
      ))
    }
    
    # Mean products
    GMP <- GSP / DFG                  # Sub-plot (genotype) mean product
    EMP_MAIN <- ESP_MAIN / DFE_MAIN   # Main plot error mean product
    EMP <- ESP / DFE                  # Sub-plot error mean product
    
    # Return for mean_products
    if (calc_type == "mean_products") {
      return(list(
        GMP = GMP,
        EMP = EMP,
        EMP_MAIN = EMP_MAIN,
        DFR = DFR,
        DFM = DFM,
        DFE_MAIN = DFE_MAIN,
        DFG = DFG,
        DFIM = DFIM,
        DFE = DFE,
        n_replications = r,
        n_main_plots = a,
        n_genotypes = b,
        design_type = "SPD"
      ))
    }
    
    # Return all (default)
    return(list(
      CF = CF,
      TSP = TSP,
      RSP = RSP,
      MSP = MSP,
      GSP = GSP,
      IMSP = IMSP,
      ESP_MAIN = ESP_MAIN,
      ESP = ESP,
      GMP = GMP,
      EMP = EMP,
      EMP_MAIN = EMP_MAIN,
      DFR = DFR,
      DFM = DFM,
      DFE_MAIN = DFE_MAIN,
      DFG = DFG,
      DFIM = DFIM,
      DFE = DFE,
      n_replications = r,
      n_main_plots = a,
      n_genotypes = b,
      design_type = "SPD"
    ))
  }
}
