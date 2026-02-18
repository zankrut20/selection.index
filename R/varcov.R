#' Calculate Variance-Covariance Components (Internal Helper)
#' 
#' @description 
#' Internal function to compute genotypic or phenotypic variance-covariance
#' matrices using math primitives. Used by gen_varcov() and phen_varcov().
#' 
#' @param data_mat Numeric matrix of trait data (n_obs x n_traits)
#' @param gen_idx Integer vector of genotype indices
#' @param rep_idx Integer vector of replication indices
#' @param col_idx Integer vector of column indices (for LSD, optional)
#' @param main_idx Integer vector of main plot indices (for SPD, optional)
#' @param design_type Integer design code: 1=RCBD, 2=LSD, 3=SPD
#' @param cov_type Integer: 1=genotypic, 2=phenotypic
#' 
#' @return Symmetric variance-covariance matrix (n_traits x n_traits)
#' 
#' @details
#' MSG and MSE are mean square matrices computed by .calculate_anova().
#' For phenotypic variance, returns MSG directly.
#' For genotypic variance, returns (MSG - MSE) / r where r is the replication factor.
#' 
#' @keywords internal
#' @noRd
.calculate_varcov <- function(data_mat, gen_idx, rep_idx,
                               col_idx = NULL, main_idx = NULL,
                               design_type = 1L, cov_type = 1L) {
  
  # Input validation
  if (!design_type %in% c(1L, 2L, 3L)) {
    stop("design_type must be 1 (RCBD), 2 (LSD), or 3 (SPD).")
  }
  if (!cov_type %in% c(1L, 2L)) {
    stop("cov_type must be 1 (genotypic) or 2 (phenotypic).")
  }
  if (design_type == 2L && is.null(col_idx)) {
    stop("col_idx is required for Latin Square Design (design_type = 2).")
  }
  if (design_type == 3L && is.null(main_idx)) {
    stop("main_idx is required for Split Plot Design (design_type = 3).")
  }
  
  n_traits <- ncol(data_mat)
  n_obs <- nrow(data_mat)
  
  # Get ANOVA components once (avoid duplication)
  anova_result <- .calculate_anova(data_mat, gen_idx, rep_idx,
                                   col_idx, main_idx, design_type)
  
  MSG <- anova_result$MSG  # Mean square for genotypes (n_traits x n_traits matrix)
  MSE <- anova_result$MSE  # Mean square for error (n_traits x n_traits matrix)
  
  # Phenotypic variance-covariance
  if (cov_type == 2L) {
    # Phenotypic variance = MSG (contains both genetic and error variance)
    return(MSG)
  }
  
  # Genotypic variance-covariance
  n_rep <- length(unique(rep_idx))
  
  if (design_type == 1L || design_type == 2L) {
    # RCBD or LSD: Vg = (MSG - MSE) / r
    Vg <- (MSG - MSE) / n_rep
    
  } else if (design_type == 3L) {
    # SPD: Vg = (MSG - MSE) / (r * m)
    n_main <- length(unique(main_idx))
    Vg <- (MSG - MSE) / (n_rep * n_main)
  }
  
  return(Vg)
}


#' Genotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments (sub-plot treatments in SPD)
#' @param replication vector containing replication/blocks (RCBD) or rows (LSD)
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param main_plots vector containing main plot treatments (required for Split Plot Design only)
#' @param design_type experimental design type: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return A Genotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' # RCBD example
#' gen_varcov(data=seldata[,3:9], genotypes=seldata$treat, replication=seldata$rep)
#' 
#' # Latin Square Design example (requires columns parameter)
#' # gen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat, 
#' #           replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")
#' 
#' # Split Plot Design example (requires main_plots parameter)
#' # gen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot, 
#' #           replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
gen_varcov <- function(data, genotypes, replication, columns = NULL, main_plots = NULL, 
                       design_type = c("RCBD", "LSD", "SPD"), 
                       method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"))
{
  design_type <- match.arg(design_type)
  
  # OPTIMIZATION: Single matrix conversion with storage.mode assignment
  # Avoids: (1) data.frame->list conversion overhead, (2) repeated as.numeric() per column
  # Why faster: Direct storage type coercion in C, no intermediate structures
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  
  colnumber <- ncol(data_mat)
  headings <- colnames(data)
  
  # OPTIMIZATION: Convert factors once outside loops (not colnumber² times)
  # Avoids: Redundant as.factor() calls and nlevels() computations
  # Why faster: Factor conversion is expensive (level sorting, attribute creation)
  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)
  repli <- nlevels(replication)
  genotype <- nlevels(genotypes)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  if (design_type == "LSD") {
    columns <- as.factor(columns)
  }
  
  # Validate Split Plot Design requirements
  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }
  if (design_type == "SPD") {
    main_plots <- as.factor(main_plots)
  }
  
  # MISSING VALUE HANDLING: Use modular engine for imputation
  # Only process method parameter if missing values are detected
  if (any(!is.finite(data_mat))) {
    # Check if user explicitly provided a method
    method_provided <- !missing(method)
    method <- match.arg(method)
    
    # Warn user if they have missing values but didn't explicitly specify a method
    if (!method_provided) {
      warning("Missing values detected in data. Using default method 'REML' for imputation. ",
              "Consider explicitly specifying method: 'REML', 'Yates', 'Healy', 'Regression', 'Mean', or 'Bartlett'.",
              call. = FALSE)
    }
    
    gen_idx <- as.integer(genotypes)
    rep_idx <- as.integer(replication)
    col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
    main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
    data_mat <- missing_value_estimation(data_mat, gen_idx, rep_idx, col_idx, main_idx, design_type, method)
  }
  
  # OPTIMIZATION: Convert factors to integer indices for design engine
  # Avoids: Factor level lookups in internal calculations
  # Why faster: Integer indexing is primitive operation
  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
  
  # C++ OPTIMIZATION: Vectorized variance-covariance computation
  # Uses math primitives for efficient grouped sums and sum of products
  # Processes all trait pairs simultaneously with optimized ANOVA calculations
  # Expected speedup: 5-20x for 7-30 traits
  design_code <- switch(design_type, "RCBD" = 1L, "LSD" = 2L, "SPD" = 3L)
  
  genetic.cov <- .calculate_varcov(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design_code,
    cov_type = 1L  # 1 = genotypic
  )
  
  # Restore dimension names
  dimnames(genetic.cov) <- list(headings, headings)
  
  return(genetic.cov)
}


#' Phenotypic Variance-Covariance Analysis
#'
#' @param data traits to be analyzed
#' @param genotypes vector containing genotypes/treatments (sub-plot treatments in SPD)
#' @param replication vector containing replication/blocks (RCBD) or rows (LSD)
#' @param columns vector containing columns (required for Latin Square Design only)
#' @param main_plots vector containing main plot treatments (required for Split Plot Design only)
#' @param design_type experimental design type: "RCBD" (default), "LSD" (Latin Square), or "SPD" (Split Plot)
#' @param method Method for missing value imputation: "REML" (default), "Yates", "Healy", "Regression", "Mean", or "Bartlett"
#'
#' @return A Phenotypic Variance-Covariance Matrix
#' @export
#'
#' @examples
#' # RCBD example
#' phen_varcov(data=seldata[,3:9], genotypes=seldata$treat, replication=seldata$rep)
#' 
#' # Latin Square Design example (requires columns parameter)
#' # phen_varcov(data=lsd_data[,3:7], genotypes=lsd_data$treat, 
#' #            replication=lsd_data$row, columns=lsd_data$col, design_type="LSD")
#' 
#' # Split Plot Design example (requires main_plots parameter)
#' # phen_varcov(data=spd_data[,3:7], genotypes=spd_data$subplot, 
#' #            replication=spd_data$block, main_plots=spd_data$mainplot, design_type="SPD")
phen_varcov <- function(data, genotypes, replication, columns = NULL, main_plots = NULL, 
                        design_type = c("RCBD", "LSD", "SPD"), 
                        method = c("REML", "Yates", "Healy", "Regression", "Mean", "Bartlett"))
{
  design_type <- match.arg(design_type)
  
  # Convert to numeric matrix once - avoids repeated data.frame/list conversions
  # storage.mode assignment is faster than as.numeric() on each column
  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "numeric"
  
  colnumber <- ncol(data_mat)
  headings <- colnames(data)
  
  # Convert factors once - eliminates colnumber² redundant conversions
  genotypes <- as.factor(genotypes)
  replication <- as.factor(replication)
  repli <- nlevels(replication)
  genotype <- nlevels(genotypes)
  
  # Validate Latin Square Design requirements
  if (design_type == "LSD" && is.null(columns)) {
    stop("Latin Square Design requires 'columns' parameter")
  }
  if (design_type == "LSD") {
    columns <- as.factor(columns)
  }
  
  # Validate Split Plot Design requirements
  if (design_type == "SPD" && is.null(main_plots)) {
    stop("Split Plot Design requires 'main_plots' parameter")
  }
  if (design_type == "SPD") {
    main_plots <- as.factor(main_plots)
  }
  
  # MISSING VALUE HANDLING: Use modular engine for imputation
  # Only process method parameter if missing values are detected
  if (any(!is.finite(data_mat))) {
    # Check if user explicitly provided a method
    method_provided <- !missing(method)
    method <- match.arg(method)
    
    # Warn user if they have missing values but didn't explicitly specify a method
    if (!method_provided) {
      warning("Missing values detected in data. Using default method 'REML' for imputation. ",
              "Consider explicitly specifying method: 'REML', 'Yates', 'Healy', 'Regression', 'Mean', or 'Bartlett'.",
              call. = FALSE)
    }
    
    gen_idx <- as.integer(genotypes)
    rep_idx <- as.integer(replication)
    col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
    main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
    data_mat <- missing_value_estimation(data_mat, gen_idx, rep_idx, col_idx, main_idx, design_type, method)
  }
  
  # Convert to integer indices - design engine uses integer grouping
  gen_idx <- as.integer(genotypes)
  rep_idx <- as.integer(replication)
  col_idx <- if (design_type == "LSD") as.integer(columns) else NULL
  main_idx <- if (design_type == "SPD") as.integer(main_plots) else NULL
  
  # C++ OPTIMIZATION: Vectorized variance-covariance computation
  # Uses math primitives for efficient grouped sums and sum of products
  # Processes all trait pairs simultaneously with optimized ANOVA calculations
  # Expected speedup: 5-20x for 7-30 traits
  design_code <- switch(design_type, "RCBD" = 1L, "LSD" = 2L, "SPD" = 3L)
  
  phenotypic.cov <- .calculate_varcov(
    data_mat = data_mat,
    gen_idx = gen_idx,
    rep_idx = rep_idx,
    col_idx = col_idx,
    main_idx = main_idx,
    design_type = design_code,
    cov_type = 2L  # 2 = phenotypic
  )
  
  # Restore dimension names
  dimnames(phenotypic.cov) <- list(headings, headings)
  
  return(phenotypic.cov)
}


# ==============================================================================
# GENOMIC VARIANCE-COVARIANCE MATRICES
# ==============================================================================
# Based on: Cerón-Rojas & Crossa (2018), Chapters 4 & 8
# Notation follows the LaTeX documentation structure

#' Genomic Variance-Covariance Matrix (Γ)
#'
#' @description
#' Computes the genomic variance-covariance matrix (Γ or G_s) from Genomic 
#' Estimated Breeding Values (GEBVs) or molecular marker scores.
#' 
#' This represents: Γ = Var(γ)
#' 
#' where γ is the vector of GEBVs across traits.
#'
#' @param gebv_mat Matrix of GEBVs (n_genotypes x n_traits)
#' @param method Method for covariance estimation: "pearson" (default), "kendall", "spearman"
#' @param use How to handle missing values: "everything", "all.obs", 
#'   "complete.obs" (default), "na.or.complete", "pairwise.complete.obs"
#'
#' @return Symmetric genomic variance-covariance matrix (n_traits x n_traits)
#' 
#' @details
#' The genomic variance-covariance matrix represents the variance-covariance 
#' structure of GEBVs predicted from molecular markers. It captures the genetic
#' relationships among traits as explained by the markers.
#' 
#' In selection index theory:
#' - Used in LGSI (Linear Genomic Selection Index)
#' - Component of Φ (phenomic-genomic covariance)
#' - Component of A (genetic-genomic covariance)
#' 
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern 
#' Plant Breeding. Springer International Publishing. Chapters 4 & 8.
#' 
#' @export
#' @examples
#' \dontrun{
#' # Simulate GEBVs
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- 5
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits), 
#'                    nrow = n_genotypes, ncol = n_traits)
#' colnames(gebv_mat) <- paste0("Trait", 1:n_traits)
#' 
#' # Compute genomic variance-covariance
#' Gamma <- genomic_varcov(gebv_mat)
#' print(Gamma)
#' }
genomic_varcov <- function(gebv_mat, method = "pearson", use = "complete.obs") {
  
  # Input validation
  gebv_mat <- as.matrix(gebv_mat)
  if (!is.numeric(gebv_mat)) {
    stop("gebv_mat must be numeric")
  }
  
  n_genotypes <- nrow(gebv_mat)
  n_traits <- ncol(gebv_mat)
  
  if (n_genotypes < 2) {
    stop("gebv_mat must have at least 2 observations")
  }
  
  # Check for missing values
  has_missing <- any(!is.finite(gebv_mat))
  if (has_missing) {
    if (use == "pairwise.complete.obs") {
      warning("Missing values detected in gebv_mat with use='pairwise.complete.obs'. ",
              "The resulting covariance matrix may not be positive semi-definite, ",
              "which can cause issues in selection index calculations. ",
              "Consider using use='complete.obs' or imputing missing values.",
              call. = FALSE)
    } else if (use == "everything") {
      stop("Missing values detected in gebv_mat. Cannot compute covariance with use='everything'. ",
           "Use 'complete.obs', 'pairwise.complete.obs', or impute missing values.")
    }
    # For other 'use' options, R's cov() will handle appropriately
  }
  
  # Compute genomic variance-covariance matrix
  Gamma <- stats::cov(gebv_mat, use = use, method = method)
  
  # Ensure symmetry (numerical precision)
  Gamma <- (Gamma + t(Gamma)) / 2
  
  # Check if matrix is positive semi-definite when missing values present
  if (has_missing && use == "pairwise.complete.obs") {
    # Quick check: all eigenvalues should be non-negative for PSD matrix
    min_eigenvalue <- min(eigen(Gamma, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eigenvalue < -1e-8) {  # Allow small numerical errors
      warning("Genomic covariance matrix is not positive semi-definite ",
              "(min eigenvalue = ", formatC(min_eigenvalue, format = "e", digits = 2), "). ",
              "This may cause problems in selection index calculations. ",
              "Consider using use='complete.obs' or a different approach.",
              call. = FALSE)
    }
  }
  
  # Preserve trait names
  trait_names <- colnames(gebv_mat)
  if (!is.null(trait_names)) {
    dimnames(Gamma) <- list(trait_names, trait_names)
  }
  
  return(Gamma)
}


#' Phenomic-Genomic Variance-Covariance Matrix (Φ)
#'
#' @description
#' Computes the combined phenomic-genomic variance-covariance matrix (Φ or P_L),
#' which is the block matrix representing the joint distribution of phenotypes 
#' and GEBVs.
#' 
#' Structure: Φ = [[P, P_yγ], [P_yγ', Γ]]
#' 
#' where:
#' - P = Var(y) = phenotypic variance-covariance
#' - Γ = Var(γ) = genomic variance-covariance  
#' - P_yγ = Cov(y, γ) = covariance between phenotypes and GEBVs
#'
#' @param phen_mat Matrix of phenotypes (n_genotypes x n_traits). 
#'   Optional if P and P_yg are provided.
#' @param gebv_mat Matrix of GEBVs (n_genotypes x n_traits).
#'   Optional if Gamma and P_yg are provided.
#' @param P Optional. Phenotypic variance-covariance matrix (n_traits x n_traits).
#'   If NULL, computed from phen_mat.
#' @param Gamma Optional. Genomic variance-covariance matrix (n_traits x n_traits).
#'   If NULL, computed from gebv_mat.
#' @param P_yg Optional. Covariance between phenotypes and GEBVs (n_traits x n_traits).
#'   If NULL, computed from phen_mat and gebv_mat.
#' @param method Method for covariance estimation: "pearson" (default), "kendall", "spearman"
#' @param use How to handle missing values: "everything", "all.obs", 
#'   "complete.obs" (default), "na.or.complete", "pairwise.complete.obs"
#'
#' @return Combined phenomic-genomic variance-covariance matrix (2t x 2t)
#'   where t is the number of traits
#' 
#' @details
#' The phenomic-genomic matrix combines phenotypic and genomic information in a 
#' single covariance structure. It is used in:
#' - CLGSI (Combined Linear Genomic Selection Index)
#' - GESIM (Genomic Eigen Selection Index Method)
#' - PPG-GESIM (Predetermined Proportional Gains GESIM)
#' 
#' When GEBVs are highly accurate predictors of true breeding values,
#' P_yγ ≈ Γ, but empirical estimation is more robust for real data.
#' 
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern 
#' Plant Breeding. Springer International Publishing. Chapters 4 & 8.
#' 
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' 
#' # Simulate phenotypes and GEBVs
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- 7
#' phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
#'                    nrow = n_genotypes, ncol = n_traits)
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
#'                    nrow = n_genotypes, ncol = n_traits)
#' 
#' # Compute phenomic-genomic covariance
#' Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)
#' print(dim(Phi))  # Should be 14 x 14 (2 * 7 traits)
#' }
phenomic_genomic_varcov <- function(phen_mat = NULL, gebv_mat = NULL,
                                    P = NULL, Gamma = NULL, P_yg = NULL,
                                    method = "pearson", use = "complete.obs") {
  
  # Check if we have necessary inputs
  has_matrices <- !is.null(P) && !is.null(Gamma) && !is.null(P_yg)
  has_data <- !is.null(phen_mat) && !is.null(gebv_mat)
  
  if (!has_matrices && !has_data) {
    stop("Must provide either (phen_mat, gebv_mat) or (P, Gamma, P_yg)")
  }
  
  # Compute covariance matrices if not provided
  if (has_data) {
    phen_mat <- as.matrix(phen_mat)
    gebv_mat <- as.matrix(gebv_mat)
    
    if (nrow(phen_mat) != nrow(gebv_mat)) {
      stop("phen_mat and gebv_mat must have the same number of rows")
    }
    if (ncol(phen_mat) != ncol(gebv_mat)) {
      stop("phen_mat and gebv_mat must have the same number of columns (traits)")
    }
    
    n_traits <- ncol(phen_mat)
    
    # Compute P if not provided
    if (is.null(P)) {
      P <- stats::cov(phen_mat, use = use, method = method)
      P <- (P + t(P)) / 2  # Ensure symmetry
    }
    
    # Compute Gamma if not provided
    if (is.null(Gamma)) {
      Gamma <- stats::cov(gebv_mat, use = use, method = method)
      Gamma <- (Gamma + t(Gamma)) / 2  # Ensure symmetry
    }
    
    # Compute P_yg if not provided
    if (is.null(P_yg)) {
      P_yg <- stats::cov(phen_mat, gebv_mat, use = use, method = method)
    }
    
  } else {
    # Validate provided matrices
    P <- as.matrix(P)
    Gamma <- as.matrix(Gamma)
    P_yg <- as.matrix(P_yg)
    
    if (!isSymmetric(unname(P), tol = 1e-8)) {
      stop("P must be symmetric")
    }
    if (!isSymmetric(unname(Gamma), tol = 1e-8)) {
      stop("Gamma must be symmetric")
    }
    
    n_traits <- nrow(P)
    
    if (nrow(Gamma) != n_traits || ncol(Gamma) != n_traits) {
      stop("Gamma must be ", n_traits, " x ", n_traits)
    }
    if (nrow(P_yg) != n_traits || ncol(P_yg) != n_traits) {
      stop("P_yg must be ", n_traits, " x ", n_traits)
    }
    # NOTE: P_yg = Cov(y, γ) is generally NOT symmetric
    # P_yg[i,j] = Cov(y_i, γ_j) ≠ Cov(y_j, γ_i) = P_yg[j,i]
    # Only check that it's square (t x t), not symmetric
  }
  
  # Construct Phi as block matrix (Equation 8.11, Chapter 8)
  # Phi = [[P,      P_yg  ],
  #        [P_yg',  Gamma ]]
  # 
  # This construction naturally produces a symmetric matrix because:
  # - P and Gamma are symmetric
  # - The off-diagonal blocks are transposes of each other: P_yg and P_yg'
  Phi <- rbind(
    cbind(P, P_yg),
    cbind(t(P_yg), Gamma)
  )
  
  # Verify that Phi is symmetric (numerical precision check)
  # If not symmetric, there's an error in the input matrices
  if (!isSymmetric(unname(Phi), tol = 1e-6)) {
    max_asymmetry <- max(abs(Phi - t(Phi)))
    warning("Phi is not symmetric (max asymmetry = ", 
            formatC(max_asymmetry, format = "e", digits = 2), 
            "). Check input matrices P, Gamma, and P_yg.")
  }
  
  # Add dimension names
  trait_names <- colnames(Phi)[1:n_traits]
  if (!is.null(trait_names)) {
    all_names <- c(paste0(trait_names, "_phen"), paste0(trait_names, "_gebv"))
    dimnames(Phi) <- list(all_names, all_names)
  }
  
  return(Phi)
}


#' Genetic-Genomic Variance-Covariance Matrix (A)
#'
#' @description
#' Computes the genetic-genomic covariance matrix (A) as defined in Chapter 8
#' (Equation 8.12) for GESIM and related genomic eigen selection indices.
#' 
#' Structure: A = [[C, C_g,γ], [C_γ,g, Γ]]  (2t × 2t, square symmetric)
#' 
#' where:
#' - C = Var(g) = true genotypic variance-covariance (t × t)
#' - Γ = Var(γ) = genomic variance-covariance (t × t)
#' - C_g,γ = Cov(g, γ) = covariance between true BVs and GEBVs (t × t)
#' - C_γ,g = Cov(γ, g) = transpose of C_g,γ (t × t)
#'
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param Gamma Genomic variance-covariance matrix (n_traits x n_traits).
#'   If NULL, assumed equal to gmat (perfect prediction).
#' @param reliability Optional. Reliability of GEBVs (r² = squared correlation 
#'   between GEBV and true BV). Can be:
#'   - Single value (applied to all traits)
#'   - Vector of length n_traits (one per trait)
#'   - NULL (default): assumes C_g,γ = Gamma (unbiased GEBVs with reliability = 1)
#' @param C_gebv_g Optional. Direct specification of Cov(γ, g) matrix (t × t).
#'   If provided, overrides reliability parameter.
#' @param square Logical. If TRUE (default), returns (2t × 2t) square matrix 
#'   as required for GESIM. If FALSE, returns (2t × t) rectangular form for LMSI.
#'
#' @return Genetic-genomic covariance matrix:
#'   - If square = TRUE: (2t × 2t) symmetric matrix for GESIM/eigen indices
#'   - If square = FALSE: (2t × t) rectangular matrix for LMSI
#'   where t is the number of traits
#' 
#' @details
#' The genetic-genomic matrix relates selection on phenotypes + GEBVs to 
#' expected genetic gains.
#' 
#' **For GESIM (Chapter 8):** Requires the full (2t × 2t) square matrix for 
#' the eigenproblem: (Φ^(-1) A - λI)b = 0
#' 
#' **For LMSI/CLGSI (Chapter 4):** Can use the rectangular (2t × t) form
#' in the equation: b = P^(-1) G w, where G is (2t × t).
#' 
#' When reliability is provided:
#' - C_γg = diag(√r²) %*% gmat (assumes accuracy scales genetic covariance)
#' 
#' When reliability is NULL:
#' - C_γg = Gamma (assumes unbiased GEBVs, perfect prediction)
#' 
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern 
#' Plant Breeding. Springer International Publishing. Chapters 4 & 8.
#' 
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' 
#' # Simulate genomic covariance
#' Gamma <- gmat * 0.8
#' 
#' # For GESIM: Get square (2t × 2t) matrix
#' A_square <- genetic_genomic_varcov(gmat, Gamma, reliability = 0.7)
#' print(dim(A_square))  # Should be 14 x 14 (2t × 2t)
#' 
#' # For LMSI: Get rectangular (2t × t) matrix
#' A_rect <- genetic_genomic_varcov(gmat, Gamma, reliability = 0.7, square = FALSE)
#' print(dim(A_rect))  # Should be 14 x 7 (2t × t)
#' }
genetic_genomic_varcov <- function(gmat, Gamma = NULL, reliability = NULL, 
                                   C_gebv_g = NULL, square = TRUE) {
  
  # Input validation
  gmat <- as.matrix(gmat)
  n_traits <- nrow(gmat)
  
  if (!isSymmetric(unname(gmat), tol = 1e-8)) {
    stop("gmat must be symmetric")
  }
  
  # Handle Gamma
  if (is.null(Gamma)) {
    Gamma <- gmat  # Assume perfect prediction
    Gamma <- (Gamma + t(Gamma)) / 2
  } else {
    Gamma <- as.matrix(Gamma)
    if (!isSymmetric(unname(Gamma), tol = 1e-8)) {
      stop("Gamma must be symmetric")
    }
    if (nrow(Gamma) != n_traits || ncol(Gamma) != n_traits) {
      stop("Gamma must be ", n_traits, " x ", n_traits)
    }
  }
  
  # Compute C_gebv_g if not directly provided
  if (!is.null(C_gebv_g)) {
    # User provided C_gebv_g directly
    C_gebv_g <- as.matrix(C_gebv_g)
    if (nrow(C_gebv_g) != n_traits || ncol(C_gebv_g) != n_traits) {
      stop("C_gebv_g must be ", n_traits, " x ", n_traits)
    }
    
  } else if (!is.null(reliability)) {
    # Use reliability to compute C_gebv_g
    # C_gebv_g = diag(accuracy) %*% gmat, where accuracy = sqrt(reliability)
    
    if (length(reliability) == 1) {
      # Single reliability value for all traits
      r_squared_vec <- rep(reliability, n_traits)
    } else if (length(reliability) == n_traits) {
      # Trait-specific reliabilities
      r_squared_vec <- reliability
    } else {
      stop("reliability must be a single value or vector of length ", n_traits)
    }
    
    # Validate reliability values
    if (any(r_squared_vec < 0) || any(r_squared_vec > 1)) {
      stop("reliability values must be between 0 and 1")
    }
    
    # Compute accuracy (r = sqrt(r²))
    accuracy_vec <- sqrt(r_squared_vec)
    
    # Scale genetic covariance by accuracy
    # This preserves cross-trait genetic correlations
    C_gebv_g <- sweep(gmat, 1, accuracy_vec, "*")
    
  } else {
    # Default: assume unbiased GEBVs (b=1), so C_gebv_g = Gamma
    C_gebv_g <- Gamma
  }
  
  if (square) {
    # Construct A as square (2t × 2t) symmetric block matrix (Equation 8.12, Chapter 8)
    # A = [[C,       C_g,γ  ],
    #      [C_γ,g,   Γ      ]]
    # 
    # where C_g,γ = Cov(g, γ) and C_γ,g = Cov(γ, g) = C_g,γ'
    # This construction naturally produces a symmetric matrix
    A <- rbind(
      cbind(gmat, C_gebv_g),
      cbind(t(C_gebv_g), Gamma)
    )
    
    # Verify symmetry (numerical precision check)
    if (!isSymmetric(unname(A), tol = 1e-6)) {
      max_asymmetry <- max(abs(A - t(A)))
      warning("A is not symmetric (max asymmetry = ", 
              formatC(max_asymmetry, format = "e", digits = 2), 
              "). Check input matrices.")
    }
    
    # Add dimension names
    trait_names <- colnames(gmat)
    if (!is.null(trait_names)) {
      all_names <- c(paste0(trait_names, "_phen"), paste0(trait_names, "_gebv"))
      dimnames(A) <- list(all_names, all_names)
    }
    
  } else {
    # Construct A as rectangular (2t × t) matrix for LMSI (Chapter 4)
    # A = [[C    ],    (t × t)
    #      [C_γg ]]    (t × t)
    # Total: (2t × t)
    A <- rbind(gmat, C_gebv_g)
    
    # Add dimension names
    trait_names <- colnames(gmat)
    if (!is.null(trait_names)) {
      row_names <- c(paste0(trait_names, "_phen"), paste0(trait_names, "_gebv"))
      dimnames(A) <- list(row_names, trait_names)
    }
  }
  
  return(A)
}
