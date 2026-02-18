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
#' @noRd
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
  
  # C++ OPTIMIZATION: Use generic math primitives
  # This allows extending to new designs without modifying C++ code
  # All design-specific logic stays in R
  
  # Ensure numeric vectors
  if (!is.numeric(trait1)) trait1 <- as.numeric(trait1)
  if (!is.numeric(trait2)) trait2 <- as.numeric(trait2)
  storage.mode(trait1) <- "numeric"
  storage.mode(trait2) <- "numeric"
  
  # Count levels once
  n_genotypes <- length(unique(genotypes))
  n_replications <- length(unique(replications))
  
  # Create data matrix for C++ functions
  data_mat <- cbind(trait1, trait2)
  
  if (design_type == "RCBD") {
    # ========== RCBD CALCULATIONS ==========
    n_obs <- n_genotypes * n_replications
    
    # Degrees of freedom
    DFG <- n_genotypes - 1
    DFR <- n_replications - 1
    DFE <- DFG * DFR
    
    # C++ PRIMITIVE: Compute grouped sums for both traits simultaneously
    gen_sums <- grouped_sums(data_mat, genotypes)
    rep_sums <- grouped_sums(data_mat, replications)
    
    # Grand totals (column sums)
    GT1 <- sum(trait1)
    GT2 <- sum(trait2)
    
    # Correction Factor
    CF <- (GT1 * GT2) / n_obs
    
    # Return early for anova_stats
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
    
    # C++ PRIMITIVE: Sum of products
    TSP <- sum(trait1 * trait2) - CF
    GSP <- sum(gen_sums[, 1] * gen_sums[, 2]) / n_replications - CF
    RSP <- sum(rep_sums[, 1] * rep_sums[, 2]) / n_genotypes - CF
    
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
    n_obs <- t * t  # In LSD, typically t×t observations
    
    # Degrees of freedom for LSD
    DFG <- t - 1           # Treatments (genotypes)
    DFR <- t - 1           # Rows
    DFC <- t - 1           # Columns
    DFE <- (t - 1) * (t - 2)  # Error: (t-1)(t-2)
    
    # C++ PRIMITIVE: Compute grouped sums for all groups
    gen_sums <- grouped_sums(data_mat, genotypes)
    row_sums <- grouped_sums(data_mat, replications)  # rows
    col_sums <- grouped_sums(data_mat, columns)
    
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
    
    # C++ PRIMITIVE: Sum of products
    TSP <- sum(trait1 * trait2) - CF
    GSP <- sum(gen_sums[, 1] * gen_sums[, 2]) / t - CF  # Genotype/Treatment SP
    RSP <- sum(row_sums[, 1] * row_sums[, 2]) / t - CF  # Row SP
    CSP <- sum(col_sums[, 1] * col_sums[, 2]) / t - CF  # Column SP
    
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
    
    n_main_plots <- length(unique(main_plots))
    n_sub_plots <- n_genotypes  # genotypes are sub-plot treatments
    n_obs <- length(trait1)
    
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
    
    # Grand totals and Correction Factor
    GT1 <- sum(trait1)
    GT2 <- sum(trait2)
    CF <- (GT1 * GT2) / n_obs
    
    # C++ PRIMITIVE: Compute grouped sums
    rep_sums <- grouped_sums(data_mat, replications)
    main_sums <- grouped_sums(data_mat, main_plots)
    gen_sums <- grouped_sums(data_mat, genotypes)
    
    # Combined factors for interactions
    rep_main_factor <- paste(replications, main_plots, sep = "_")
    main_sub_factor <- paste(main_plots, genotypes, sep = "_")
    
    rep_main_sums <- grouped_sums(data_mat, rep_main_factor)
    main_sub_sums <- grouped_sums(data_mat, main_sub_factor)
    
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
    
    # C++ PRIMITIVE: Sum of products
    TSP <- sum(trait1 * trait2) - CF
    RSP <- sum(rep_sums[, 1] * rep_sums[, 2]) / (a * b) - CF
    MSP <- sum(main_sums[, 1] * main_sums[, 2]) / (r * b) - CF
    GSP <- sum(gen_sums[, 1] * gen_sums[, 2]) / (r * a) - CF
    RMSP <- sum(rep_main_sums[, 1] * rep_main_sums[, 2]) / b - CF
    
    # Main plot error sum of products
    ESP_MAIN <- RMSP - RSP - MSP
    
    # Main × Sub interaction sum of products
    IMSP <- sum(main_sub_sums[, 1] * main_sub_sums[, 2]) / r - CF - MSP - GSP
    
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


#' Design Statistics API - Single Engine for ANOVA Computations
#'
#' @description
#' Unified API for experimental design statistics and ANOVA computations.
#' Replaces ad-hoc ANOVA calculations throughout the package, providing
#' a single source of truth for correction factors, sums of products,
#' degrees of freedom, and mean squares.
#'
#' This function is a wrapper around design_stats() that computes
#' multivariate ANOVA statistics (mean square matrices) for all trait
#' pairs simultaneously.
#'
#' @param data_mat Numeric matrix of trait data (n_obs x n_traits)
#' @param gen_idx Integer vector of genotype indices (sub-plot treatments in SPD)
#' @param rep_idx Integer vector of replication/block indices (RCBD) or row indices (LSD)
#' @param col_idx Integer vector of column indices (for LSD, optional)
#' @param main_idx Integer vector of main plot indices (for SPD, optional)
#' @param design_type Integer design code: 1=RCBD, 2=LSD, 3=SPD
#'
#' @return List with components compatible with legacy .calculate_anova():
#'   \item{GMS}{Genotype mean squares vector (diagonal of MSG)}
#'   \item{EMS}{Error mean squares vector (diagonal of MSE)}
#'   \item{EMS_MAIN}{Main plot error mean squares vector (SPD only, diagonal of MSG_MAIN)}
#'   \item{DFG}{Degrees of freedom for genotypes/sub-plots}
#'   \item{DFE}{Degrees of freedom for error (sub-plot error for SPD)}
#'   \item{DFE_MAIN}{Degrees of freedom for main plot error (SPD only)}
#'   \item{n_rep}{Number of replications}
#'   \item{n_gen}{Number of genotypes/sub-plot treatments}
#'   \item{n_main}{Number of main plot treatments (SPD only)}
#'   \item{MSG}{Genotype mean square matrix (n_traits x n_traits)}
#'   \item{MSE}{Error mean square matrix (n_traits x n_traits)}
#'
#' @details
#' This function centralizes ANOVA computation logic, eliminating code
#' duplication across varcov.R, mean_performance.R, and other modules.
#' 
#' **Design-specific formulas:**
#' 
#' RCBD:
#' - MSG = GSP / (g - 1)
#' - MSE = ESP / ((g - 1) * (r - 1))
#' 
#' LSD:
#' - MSG = GSP / (t - 1)
#' - MSE = ESP / ((t - 1) * (t - 2))
#' 
#' SPD:
#' - MSG = GSP / (b - 1)  [sub-plot treatments]
#' - MSE = ESP / (a * (b - 1) * r)  [sub-plot error]
#' - MSE_MAIN = ESP_MAIN / ((a - 1) * (r - 1))  [main plot error]
#'
#' @keywords internal
#' @noRd
design_stats_api <- function(data_mat, gen_idx, rep_idx,
                              col_idx = NULL, main_idx = NULL,
                              design_type = 1L) {
  
  # Convert design_type integer to character
  design_char <- switch(as.character(design_type),
                        "1" = "RCBD",
                        "2" = "LSD", 
                        "3" = "SPD",
                        stop("design_type must be 1 (RCBD), 2 (LSD), or 3 (SPD)"))
  
  n_traits <- ncol(data_mat)
  
  # Initialize mean square matrices
  MSG <- matrix(0, n_traits, n_traits)
  MSE <- matrix(0, n_traits, n_traits)
  MSE_MAIN <- if (design_type == 3L) matrix(0, n_traits, n_traits) else NULL
  
  # Compute all pairwise trait statistics
  for (i in seq_len(n_traits)) {
    for (j in i:n_traits) {
      # Call design_stats for this trait pair
      stats <- design_stats(
        trait1 = data_mat[, i],
        trait2 = data_mat[, j],
        genotypes = gen_idx,
        replications = rep_idx,
        columns = col_idx,
        main_plots = main_idx,
        design_type = design_char,
        calc_type = "all"
      )
      
      # Extract mean products and populate matrices
      MSG[i, j] <- MSG[j, i] <- stats$GMP
      MSE[i, j] <- MSE[j, i] <- stats$EMP
      
      # For SPD, also populate main plot error matrix
      if (design_type == 3L && !is.null(stats$EMP_MAIN)) {
        MSE_MAIN[i, j] <- MSE_MAIN[j, i] <- stats$EMP_MAIN
      }
    }
  }
  
  # Get degrees of freedom and counts from last iteration
  # (they are the same for all trait pairs)
  final_stats <- design_stats(
    trait1 = data_mat[, 1],
    trait2 = data_mat[, 1],
    genotypes = gen_idx,
    replications = rep_idx,
    columns = col_idx,
    main_plots = main_idx,
    design_type = design_char,
    calc_type = "anova_stats"
  )
  
  # Extract vectors for diagonal elements (backward compatibility)
  GMS_vec <- diag(MSG)
  EMS_vec <- diag(MSE)
  EMS_MAIN_vec <- if (design_type == 3L) diag(MSE_MAIN) else rep(NA_real_, n_traits)
  
  # Return in same format as .calculate_anova()
  list(
    GMS = GMS_vec,           # Vector for mean_performance
    EMS = EMS_vec,           # Vector for mean_performance 
    EMS_MAIN = EMS_MAIN_vec, # Vector for SPD
    DFG = final_stats$DFG,   # Degrees of freedom genotype
    DFE = final_stats$DFE,   # Degrees of freedom error
    DFE_MAIN = if (design_type == 3L) final_stats$DFE_MAIN else NA_integer_,
    n_rep = final_stats$n_replications,
    n_gen = final_stats$n_genotypes,
    n_main = if (design_type == 3L) final_stats$n_main_plots else NA_integer_,
    MSG = MSG,               # Matrix for variance-covariance
    MSE = MSE                # Matrix for variance-covariance
  )
}
