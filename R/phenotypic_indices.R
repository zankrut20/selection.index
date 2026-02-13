#' Phenotypic Selection Indices (Chapter 2)
#' @name phenotypic_indices
#'
#' @description
#' Implements phenotypic selection index methods from Chapter 2:
#' The Linear Phenotypic Selection Index.
#'
#' Methods included:
#' - Smith-Hazel Index (LPSI) - The optimal unrestricted phenotypic index
#' - Base Index - Simple index using economic weights directly
#' - Combinatorial Index Builder - Builds indices for all trait combinations
#'
#' All implementations use C++ primitives for mathematical operations.
#'
#' @references
#' Smith, H. F. (1936). A discriminant function for plant selection.
#' Annals of Eugenics, 7(3), 240-250.
#'
#' Hazel, L. N. (1943). The genetic basis for constructing selection indexes.
#' Genetics, 28(6), 476.
#'
#' Williams, J. S. (1962). The evaluation of a selection index.
#' Biometrics, 18(3), 375-393.
#'
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 2.
#'
#' @keywords internal
#' @importFrom stats cov sd cor setNames
#' @importFrom utils combn
NULL

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Solve symmetric system for multiple right-hand sides
#' @keywords internal
.solve_sym_multi <- function(A, B) {
  B <- as.matrix(B)
  n_col <- ncol(B)
  out <- matrix(0, nrow = nrow(A), ncol = n_col)
  for (j in seq_len(n_col)) {
    out[, j] <- cpp_symmetric_solve(A, B[, j])
  }
  out
}

#' Compute selection index metrics using C++ primitives
#' @keywords internal
.index_metrics <- function(b, P, G, w = NULL, const_factor = 2.063, GAY = NULL) {
  b <- as.numeric(b)
  
  # Variance of index: σ²_I = b'Pb
  bPb <- cpp_quadratic_form_sym(b, P)
  
  # Genetic variance of index: b'Gb
  bGb <- cpp_quadratic_form_sym(b, G)
  
  # Standard deviation of index
  sigma_I <- if (bPb > 0) sqrt(bPb) else NA_real_
  
  # Total genetic advance: R_H = i * σ_I
  delta_g_scalar <- if (!is.na(sigma_I)) const_factor * sigma_I else NA_real_
  
  # Expected genetic gain per trait: ΔG = (i/σ_I) * Gb
  delta_g_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    const_factor * (G %*% b) / sigma_I
  } else {
    rep(NA_real_, nrow(G))
  }
  
  # Heritability of index: h²_I = b'Gb / b'Pb
  hI2 <- if (!is.na(bPb) && bPb > 0) bGb / bPb else NA_real_
  
  # Accuracy of index: r_HI = sqrt(b'Gb / b'Pb)
  rHI <- if (!is.na(hI2) && hI2 >= 0) sqrt(hI2) else NA_real_
  
  # Overall genetic advance (if weights provided)
  GA <- NA_real_
  PRE <- NA_real_
  if (!is.null(w)) {
    bGw <- cpp_quadratic_form(b, G, w)
    GA <- if (!is.na(sigma_I) && sigma_I > 0) const_factor * bGw / sigma_I else NA_real_
    PRE_constant <- if (is.null(GAY)) 100 else 100 / GAY
    PRE <- if (!is.na(GA)) GA * PRE_constant else NA_real_
  }
  
  list(
    bPb = bPb,
    bGb = bGb,
    sigma_I = sigma_I,
    Delta_G = delta_g_scalar,
    Delta_G_vec = as.vector(delta_g_vec),
    hI2 = hI2,
    rHI = rHI,
    GA = GA,
    PRE = PRE
  )
}

# ==============================================================================
# SMITH-HAZEL LINEAR PHENOTYPIC SELECTION INDEX (LPSI)
# ==============================================================================

#' Smith-Hazel Linear Phenotypic Selection Index
#'
#' @description
#' Implements the optimal Smith-Hazel selection index which maximizes
#' the correlation between the index I = b'y and the breeding objective H = w'g.
#'
#' This is the foundational selection index method from Chapter 2.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param wmat Economic weights matrix (n_traits x k), or vector
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param selection_intensity Selection intensity constant (default: 2.063 for 10% selection)
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients and metrics
#'     \item \code{b} - Vector of Smith-Hazel index coefficients
#'     \item \code{w} - Named vector of economic weights
#'     \item \code{Delta_G} - Named vector of expected genetic gains per trait
#'     \item \code{sigma_I} - Standard deviation of the index
#'     \item \code{GA} - Total genetic advance
#'     \item \code{PRE} - Percent relative efficiency
#'     \item \code{hI2} - Heritability of the index
#'     \item \code{rHI} - Accuracy (correlation with breeding objective)
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 2):}
#'
#' Index coefficients: \eqn{b = P^{-1}Gw}
#'
#' Where:
#' - P = Phenotypic variance-covariance matrix
#' - G = Genotypic variance-covariance matrix
#' - w = Economic weights
#'
#' Key metrics:
#' - Variance of index: \eqn{\sigma^2_I = b'Pb}
#' - Total genetic advance: \eqn{R_H = i\sqrt{b'Pb}}
#' - Expected gains per trait: \eqn{\Delta G = (i/\sigma_I)Gb}
#' - Heritability of index: \eqn{h^2_I = b'Gb / b'Pb}
#' - Accuracy: \eqn{r_{HI} = \sqrt{b'Gb / b'Pb}}
#'
#' @export
#' @examples
#' \dontrun{
#' # Calculate variance-covariance matrices
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' 
#' # Define economic weights
#' weights <- c(10, 8, 6, 4, 2, 1, 1)
#' 
#' # Build Smith-Hazel index
#' result <- smith_hazel(pmat, gmat, weights)
#' print(result)
#' summary(result)
#' }
smith_hazel <- function(pmat, gmat, wmat,
                        wcol = 1,
                        selection_intensity = 2.063,
                        GAY = NULL) {
  
  # ============================================================================
  # STEP 1: Input Validation and Preparation
  # ============================================================================
  
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  
  n_traits <- nrow(pmat)
  
  if (nrow(pmat) != ncol(pmat) || nrow(gmat) != ncol(gmat)) {
    stop("pmat and gmat must be square matrices")
  }
  
  if (nrow(pmat) != nrow(gmat)) {
    stop("pmat and gmat must have the same dimensions")
  }
  
  # Handle wmat: convert vector to matrix if needed
  if (is.vector(wmat)) {
    wmat <- matrix(wmat, ncol = 1)
  } else {
    wmat <- as.matrix(wmat)
  }
  
  if (nrow(wmat) != n_traits) {
    stop("Number of rows in wmat must equal number of traits (", n_traits, ")")
  }
  
  if (wcol < 1 || wcol > ncol(wmat)) {
    stop("wcol must be between 1 and ", ncol(wmat))
  }
  
  # Extract economic weights from specified column
  w <- as.numeric(wmat[, wcol])
  
  if (any(!is.finite(w))) {
    stop("Economic weights must be finite")
  }
  
  # ============================================================================
  # STEP 2: Calculate Smith-Hazel Index Coefficients: b = P^{-1}Gw
  # ============================================================================
  
  # Compute Gw using matrix multiplication
  Gw <- gmat %*% w
  
  # Solve P * b = Gw using symmetric solver
  b <- cpp_symmetric_solve(pmat, Gw)
  
  if (any(!is.finite(b))) {
    stop("Failed to compute index coefficients. Check matrix conditioning.")
  }
  
  # ============================================================================
  # STEP 3: Calculate Index Metrics
  # ============================================================================
  
  metrics <- .index_metrics(
    b = b,
    P = pmat,
    G = gmat,
    w = w,
    const_factor = selection_intensity,
    GAY = GAY
  )
  
  # ============================================================================
  # STEP 4: Build Summary Data Frame
  # ============================================================================
  
  b_vec <- round(b, 4)
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))
  
  summary_df <- data.frame(
    b_df,
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    Delta_G = round(metrics$Delta_G, 4),
    hI2 = round(metrics$hI2, 4),
    rHI = round(metrics$rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # ============================================================================
  # STEP 5: Build Return Object
  # ============================================================================
  
  result <- list(
    b = as.numeric(b_vec),
    w = setNames(w, colnames(pmat)),
    Delta_G = setNames(metrics$Delta_G_vec, colnames(pmat)),
    sigma_I = metrics$sigma_I,
    GA = metrics$GA,
    PRE = metrics$PRE,
    hI2 = metrics$hI2,
    rHI = metrics$rHI,
    selection_intensity = selection_intensity,
    summary = summary_df
  )
  
  class(result) <- c("smith_hazel", "lpsi", "selection_index", "list")
  
  return(result)
}

# ==============================================================================
# BASE INDEX (Williams, 1962)
# ==============================================================================

#' Base Index (Williams, 1962)
#'
#' @description
#' Implements the Base Index where coefficients are set equal to economic weights.
#' This is a simple, non-optimized approach that serves as a baseline comparison.
#'
#' Unlike the Smith-Hazel index which requires matrix inversion, the Base Index
#' is computationally trivial and robust when covariance estimates are unreliable.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param wmat Economic weights matrix (n_traits x k), or vector
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param selection_intensity Selection intensity constant (default: 2.063)
#' @param compare_to_lpsi Logical. If TRUE, compares Base Index efficiency to optimal LPSI (default: TRUE)
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients and metrics
#'     \item \code{b} - Vector of Base Index coefficients (equal to w)
#'     \item \code{w} - Named vector of economic weights
#'     \item \code{Delta_G} - Named vector of expected genetic gains per trait
#'     \item \code{lpsi_comparison} - Optional comparison with Smith-Hazel LPSI
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' Index coefficients: \eqn{b = w}
#'
#' The Base Index is appropriate when:
#' - Covariance estimates are unreliable
#' - Computational simplicity is required
#' - A baseline for comparison is needed
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' weights <- c(10, 8, 6, 4, 2, 1, 1)
#' 
#' result <- base_index(pmat, gmat, weights, compare_to_lpsi = TRUE)
#' print(result)
#' }
base_index <- function(pmat, gmat, wmat,
                       wcol = 1,
                       selection_intensity = 2.063,
                       compare_to_lpsi = TRUE,
                       GAY = NULL) {
  
  # ============================================================================
  # STEP 1: Input Validation and Preparation
  # ============================================================================
  
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  
  n_traits <- nrow(pmat)
  
  if (nrow(pmat) != ncol(pmat) || nrow(gmat) != ncol(gmat)) {
    stop("pmat and gmat must be square matrices")
  }
  
  if (nrow(pmat) != nrow(gmat)) {
    stop("pmat and gmat must have the same dimensions")
  }
  
  # Handle wmat: convert vector to matrix if needed
  if (is.vector(wmat)) {
    wmat <- matrix(wmat, ncol = 1)
  } else {
    wmat <- as.matrix(wmat)
  }
  
  if (nrow(wmat) != n_traits) {
    stop("Number of rows in wmat must equal number of traits (", n_traits, ")")
  }
  
  if (wcol < 1 || wcol > ncol(wmat)) {
    stop("wcol must be between 1 and ", ncol(wmat))
  }
  
  # Extract economic weights from specified column
  w <- as.numeric(wmat[, wcol])
  
  if (any(!is.finite(w))) {
    stop("Economic weights must be finite")
  }
  
  # ============================================================================
  # STEP 2: Base Index Calculation (b = w)
  # ============================================================================
  
  # In Base Index, coefficients are simply the economic weights
  b <- w
  
  # ============================================================================
  # STEP 3: Calculate Index Metrics
  # ============================================================================
  
  metrics <- .index_metrics(
    b = b,
    P = pmat,
    G = gmat,
    w = w,
    const_factor = selection_intensity,
    GAY = GAY
  )
  
  # ============================================================================
  # STEP 4: Optional Comparison with LPSI
  # ============================================================================
  
  lpsi_comparison <- NULL
  if (compare_to_lpsi) {
    # Calculate optimal LPSI coefficients: b_lpsi = P^{-1} G w
    tryCatch({
      P_inv_G <- .solve_sym_multi(pmat, gmat)
      b_lpsi <- P_inv_G %*% w
      
      metrics_lpsi <- .index_metrics(
        b = b_lpsi,
        P = pmat,
        G = gmat,
        w = w,
        const_factor = selection_intensity,
        GAY = GAY
      )
      
      lpsi_comparison <- list(
        b_lpsi = as.numeric(b_lpsi),
        GA_lpsi = metrics_lpsi$GA,
        PRE_lpsi = metrics_lpsi$PRE,
        hI2_lpsi = metrics_lpsi$hI2,
        rHI_lpsi = metrics_lpsi$rHI,
        Delta_G_lpsi = setNames(metrics_lpsi$Delta_G_vec, colnames(pmat)),
        efficiency_ratio = if (!is.na(metrics$GA) && !is.na(metrics_lpsi$GA) && metrics_lpsi$GA > 0) {
          metrics$GA / metrics_lpsi$GA
        } else {
          NA_real_
        }
      )
    }, error = function(e) {
      warning("Could not calculate LPSI comparison: ", e$message, call. = FALSE)
      lpsi_comparison <- NULL
    })
  }
  
  # ============================================================================
  # STEP 5: Build Summary Data Frame
  # ============================================================================
  
  b_vec <- round(b, 4)
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))
  
  summary_df <- data.frame(
    b_df,
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    Delta_G = round(metrics$Delta_G, 4),
    hI2 = round(metrics$hI2, 4),
    rHI = round(metrics$rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # ============================================================================
  # STEP 6: Build Return Object
  # ============================================================================
  
  result <- list(
    b = b_vec,
    w = setNames(w, colnames(pmat)),
    Delta_G = setNames(metrics$Delta_G_vec, colnames(pmat)),
    sigma_I = metrics$sigma_I,
    GA = metrics$GA,
    PRE = metrics$PRE,
    hI2 = metrics$hI2,
    rHI = metrics$rHI,
    selection_intensity = selection_intensity,
    summary = summary_df,
    lpsi_comparison = lpsi_comparison
  )
  
  class(result) <- c("base_index", "selection_index", "list")
  
  return(result)
}

# ==============================================================================
# COMBINATORIAL LPSI (All Trait Combinations)
# ==============================================================================

#' Combinatorial Linear Phenotypic Selection Index
#'
#' @description
#' Build all possible Smith-Hazel selection indices from trait combinations,
#' with optional exclusion of specific traits.
#'
#' This function systematically evaluates indices for all combinations of
#' ncomb traits, which is useful for identifying the most efficient subset
#' of traits for selection.
#'
#' @param ncomb Number of traits per combination
#' @param pmat Phenotypic variance-covariance matrix
#' @param gmat Genotypic variance-covariance matrix
#' @param wmat Weight matrix
#' @param wcol Weight column number if more than one weight set (default: 1)
#' @param GAY Genetic advance of comparative trait (optional)
#' @param excluding_trait Optional. Traits to exclude from combinations. Can be:
#'   (1) numeric vector of trait indices (e.g., c(1, 3)),
#'   (2) character vector of trait names (e.g., c("sypp", "dtf")),
#'   (3) data frame/matrix columns with trait data (trait names extracted from column names).
#'   When specified, only combinations that do NOT contain any of these traits are returned.
#'
#' @return Data frame of all possible selection indices with metrics (GA, PRE, Delta_G, rHI, hI2)
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' wmat <- weight_mat(weight)
#' 
#' # Build all 3-trait indices
#' result <- lpsi(ncomb = 3, pmat = pmat, gmat = gmat, wmat = wmat, wcol = 1)
#' 
#' # Exclude specific traits
#' result <- lpsi(ncomb = 3, pmat = pmat, gmat = gmat, wmat = wmat, 
#'                excluding_trait = c(1, 3))
#' }
lpsi <- function(ncomb, pmat, gmat, wmat, wcol = 1, GAY, excluding_trait = NULL){
  # Convert matrices once
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  wmat <- as.matrix(wmat)
  
  # Process excluding_trait parameter
  exclude_indices <- NULL
  if (!is.null(excluding_trait)) {
    # Case 1: numeric vector of trait indices
    if (is.numeric(excluding_trait)) {
      exclude_indices <- unique(as.integer(excluding_trait))
    }
    # Case 2: character vector of trait names
    else if (is.character(excluding_trait)) {
      trait_names <- colnames(pmat)
      if (is.null(trait_names)) {
        stop("pmat must have column names to use character trait names in excluding_trait")
      }
      exclude_indices <- which(trait_names %in% excluding_trait)
      if (length(exclude_indices) == 0) {
        warning("None of the specified trait names found in pmat column names")
      }
    }
    # Case 3: data frame or matrix (extract column names)
    else if (is.data.frame(excluding_trait) || is.matrix(excluding_trait)) {
      exclude_names <- colnames(excluding_trait)
      if (is.null(exclude_names)) {
        stop("excluding_trait data must have column names")
      }
      trait_names <- colnames(pmat)
      if (is.null(trait_names)) {
        stop("pmat must have column names to match with excluding_trait column names")
      }
      exclude_indices <- which(trait_names %in% exclude_names)
      if (length(exclude_indices) == 0) {
        warning("None of the column names from excluding_trait found in pmat")
      }
    }
    else {
      stop("excluding_trait must be a numeric vector, character vector, data frame, or matrix")
    }
  }
  
  # Generate combinations (keep as columns)
  ncolmn <- ncol(pmat)
  comb_all <- combn(ncolmn, ncomb)
  
  # Filter combinations if excluding_trait is specified
  if (!is.null(exclude_indices) && length(exclude_indices) > 0) {
    # Vectorized filtering: keep combinations that don't contain any excluded traits
    keep_mask <- colSums(matrix(comb_all %in% exclude_indices, nrow = nrow(comb_all))) == 0
    comb <- comb_all[, keep_mask, drop = FALSE]
    
    # Handle edge case: no combinations after filtering
    if (ncol(comb) == 0) {
      return(data.frame(ID = character(0), GA = numeric(0), 
                        PRE = numeric(0), Delta_G = numeric(0),
                        rHI = numeric(0), hI2 = numeric(0), 
                        Rank = numeric(0)))
    }
  } else {
    comb <- comb_all
  }
  
  ncomb_total <- ncol(comb)
  
  # Pre-compute constants
  const_factor <- 2.063
  PRE_constant <- if(missing(GAY)) 100 else 100 / GAY
  
  # CRITICAL: Compute FULL breeding objective ONCE (H = w'g for ALL traits)
  # The goal is to find which trait subset best predicts this FIXED H
  w_full <- cpp_extract_vector(wmat, seq_len(ncolmn), wcol - 1L)
  Gw_full <- gmat %*% w_full  # Cov(g_i, H) for all traits
  wGw_full <- cpp_quadratic_form_sym(w_full, gmat)  # Var(H) - constant denominator
  
  # Pre-allocate result storage
  IDs <- character(ncomb_total)
  b_list <- vector("list", ncomb_total)
  GAs <- numeric(ncomb_total)
  PREs <- numeric(ncomb_total)
  Delta_Gs <- numeric(ncomb_total)
  rHIs <- numeric(ncomb_total)
  hI2s <- numeric(ncomb_total)
  
  # Process each combination using math primitives
  for (j in seq_len(ncomb_total)) {
    # Get trait indices for this combination (1-indexed)
    idx <- comb[, j]
    IDs[j] <- paste(idx, collapse = ", ")
    
    # Extract submatrices for this combination
    P_sub <- cpp_extract_submatrix(pmat, idx)
    G_sub <- cpp_extract_submatrix(gmat, idx)
    
    # CRITICAL FIX: Extract Gw for subset from pre-computed FULL Gw
    # This ensures all combinations optimize for the SAME breeding objective H
    Gw_sub <- Gw_full[idx, , drop = FALSE]
    
    # Calculate selection index coefficients: b = P_sub^(-1) * Gw_sub
    # where Gw_sub = Cov(y_subset, H_full)
    b <- cpp_symmetric_solve(P_sub, Gw_sub)
    
    # Calculate quadratic forms needed for all metrics
    bPb <- cpp_quadratic_form_sym(b, P_sub)  # b'Pb (variance of index)
    bGb <- cpp_quadratic_form_sym(b, G_sub)  # b'Gb (genetic variance of index)
    
    # CRITICAL FIX: Correlation with FULL breeding objective
    # Cov(I, H) = b' * Cov(y_subset, H_full) = b' * Gw_sub
    bGw_full <- as.numeric(crossprod(b, Gw_sub))  # Cov(I, H_full)
    
    # Genetic Advance: GA = i * Cov(I, H) / sqrt(Var(I))
    sigma_I <- sqrt(bPb)
    GA <- const_factor * bGw_full / sigma_I
    
    # Percent relative efficiency
    PRE <- GA * PRE_constant
    
    # Selection response: R = i * sigma_I
    Delta_G <- const_factor * sigma_I
    
    # Index Heritability: hI² = b'Gb / b'Pb
    hI2 <- if (bPb > 0) bGb / bPb else 0
    
    # Accuracy: r_HI = Cov(I, H) / (sqrt(Var(I)) * sqrt(Var(H)))
    # Using FULL Var(H) for all combinations
    rHI <- if (bPb > 0 && wGw_full > 0) {
      abs(bGw_full) / (sigma_I * sqrt(wGw_full))
    } else {
      0
    }
    
    # Store results (rounded to 4 decimals)
    b_list[[j]] <- round(as.vector(b), 4)
    GAs[j] <- round(GA, 4)
    PREs[j] <- round(PRE, 4)
    Delta_Gs[j] <- round(Delta_G, 4)
    rHIs[j] <- round(rHI, 4)
    hI2s[j] <- round(hI2, 4)
  }
  
  # Convert b_list to matrix (pad with NA for shorter vectors)
  max_b_cols <- max(sapply(b_list, length))
  b_matrix <- matrix(NA_real_, nrow = ncomb_total, ncol = max_b_cols)
  colnames(b_matrix) <- paste0("b.", seq_len(max_b_cols))
  
  for (j in seq_len(ncomb_total)) {
    b_len <- length(b_list[[j]])
    b_matrix[j, 1:b_len] <- b_list[[j]]
  }
  
  # Construct result data frame
  df <- data.frame(
    ID = IDs,
    b_matrix,
    GA = GAs,
    PRE = PREs,
    Delta_G = Delta_Gs,
    rHI = rHIs,
    hI2 = hI2s,
    Rank = rank(-PREs, ties.method = "min"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  return(df)
}

# ==============================================================================
# PRINT AND SUMMARY METHODS
# ==============================================================================

#' Print method for Smith-Hazel Index
#'
#' @param x Object of class 'smith_hazel'
#' @param ... Additional arguments (unused)
#' @export
print.smith_hazel <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("SMITH-HAZEL LINEAR PHENOTYPIC SELECTION INDEX\n")
  cat("==============================================================\n\n")
  
  cat("Selection intensity (i):", x$selection_intensity, "\n\n")
  
  cat("Index Coefficients: b = P^{-1}Gw\n\n")
  
  # Handle missing trait names
  trait_names <- names(x$w)
  if (is.null(trait_names) || length(trait_names) == 0) {
    trait_names <- paste0("Trait_", seq_along(x$w))
  }
  
  coef_df <- data.frame(
    Trait = trait_names,
    Economic_Weight = round(as.numeric(x$w), 4),
    Coefficient = round(as.numeric(x$b), 4),
    stringsAsFactors = FALSE
  )
  print(coef_df)
  
  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("Genetic Advance (GA):         %.4f\n", x$GA))
  cat(sprintf("Index Heritability (hI2):     %.4f\n", x$hI2))
  cat(sprintf("Accuracy (r_HI):              %.4f\n", x$rHI))
  cat(sprintf("Index Std Dev (sigma_I):      %.4f\n", x$sigma_I))
  
  if (!is.na(x$PRE)) {
    cat(sprintf("Relative Efficiency (PRE):    %.2f%%\n", x$PRE))
  }
  
  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC RESPONSE PER TRAIT\n")
  cat("-------------------------------------------------------------\n")
  
  response_df <- data.frame(
    Trait = trait_names,
    Delta_G = round(as.numeric(x$Delta_G), 4),
    stringsAsFactors = FALSE
  )
  print(response_df)
  
  cat("\n")
  invisible(x)
}

#' Summary method for Smith-Hazel Index
#'
#' @param object Object of class 'smith_hazel'
#' @param ... Additional arguments passed to print
#' @export
summary.smith_hazel <- function(object, ...) {
  # Print standard output first
  print(object, ...)
  
  # Add additional summary statistics
  cat("\n")
  cat("==============================================================\n")
  cat("ADDITIONAL STATISTICS\n")
  cat("==============================================================\n\n")
  
  cat("Economic Weights:\n")
  cat(sprintf("  Mean:         %.4f\n", mean(object$w)))
  cat(sprintf("  SD:           %.4f\n", sd(object$w)))
  cat(sprintf("  Range:        [%.4f, %.4f]\n", min(object$w), max(object$w)))
  
  cat("\nExpected Genetic Gains:\n")
  cat(sprintf("  Mean:         %.4f\n", mean(object$Delta_G)))
  cat(sprintf("  SD:           %.4f\n", sd(object$Delta_G)))
  cat(sprintf("  Range:        [%.4f, %.4f]\n", min(object$Delta_G), max(object$Delta_G)))
  
  cat("\nIndex Coefficients:\n")
  cat(sprintf("  Mean:         %.4f\n", mean(object$b)))
  cat(sprintf("  SD:           %.4f\n", sd(object$b)))
  cat(sprintf("  Range:        [%.4f, %.4f]\n", min(object$b), max(object$b)))
  
  cat("\n")
  invisible(object)
}

#' Print method for Base Index
#'
#' @param x Object of class 'base_index'
#' @param ... Additional arguments (unused)
#' @export
print.base_index <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("BASE INDEX (Williams, 1962)\n")
  cat("==============================================================\n\n")
  
  cat("Selection intensity (i):", x$selection_intensity, "\n\n")
  
  cat("Index Coefficients (b = w):\n")
  cat("The Base Index sets coefficients equal to economic weights.\n\n")
  
  # Handle missing trait names
  trait_names <- names(x$w)
  if (is.null(trait_names) || length(trait_names) == 0) {
    trait_names <- paste0("Trait_", seq_along(x$w))
  }
  
  coef_df <- data.frame(
    Trait = trait_names,
    Economic_Weight = round(as.numeric(x$w), 4),
    Coefficient = round(as.numeric(x$b), 4),
    stringsAsFactors = FALSE
  )
  print(coef_df)
  
  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("Genetic Advance (GA):         %.4f\n", x$GA))
  cat(sprintf("Index Heritability (hI2):     %.4f\n", x$hI2))
  cat(sprintf("Correlation (H, I):           %.4f\n", x$rHI))
  cat(sprintf("Index Std Dev (sigma_I):      %.4f\n", x$sigma_I))
  
  if (!is.na(x$PRE)) {
    cat(sprintf("Relative Efficiency (PRE):    %.2f%%\n", x$PRE))
  }
  
  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC RESPONSE\n")
  cat("-------------------------------------------------------------\n")
  
  response_df <- data.frame(
    Trait = trait_names,
    Delta_G = round(as.numeric(x$Delta_G), 4),
    stringsAsFactors = FALSE
  )
  print(response_df)
  
  if (!is.null(x$lpsi_comparison)) {
    cat("\n\n")
    cat("-------------------------------------------------------------\n")
    cat("COMPARISON WITH OPTIMAL LPSI\n")
    cat("-------------------------------------------------------------\n")
    cat(sprintf("Base Index GA:        %.4f\n", x$GA))
    cat(sprintf("Optimal LPSI GA:      %.4f\n", x$lpsi_comparison$GA_lpsi))
    cat(sprintf("Efficiency Ratio:     %.2f%% of LPSI\n",
                x$lpsi_comparison$efficiency_ratio * 100))
    
    if (x$lpsi_comparison$efficiency_ratio < 0.9) {
      cat("\n[!] Base Index achieves <90% of LPSI efficiency.\n")
      cat("  Consider using optimized LPSI if covariance estimates are reliable.\n")
    } else if (x$lpsi_comparison$efficiency_ratio >= 0.95) {
      cat("\n[OK] Base Index performs well (>=95% of LPSI efficiency).\n")
      cat("  The simple Base Index is adequate for this scenario.\n")
    }
  }
  
  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("INTERPRETATION\n")
  cat("-------------------------------------------------------------\n")
  cat("The Base Index (b = w) is a simple, unoptimized approach that:\n")
  cat("  - Does not require matrix inversion\n")
  cat("  - Is robust when covariance estimates are unreliable\n")
  cat("  - Serves as a baseline for comparing optimized indices\n")
  cat("  - May be less efficient than LPSI but more stable\n")
  
  cat("\n")
  invisible(x)
}

#' Summary method for Base Index
#'
#' @param object Object of class 'base_index'
#' @param ... Additional arguments passed to print
#' @export
summary.base_index <- function(object, ...) {
  # Print standard output first
  print(object, ...)
  
  # Add additional summary statistics
  cat("\n")
  cat("==============================================================\n")
  cat("ADDITIONAL DETAILS\n")
  cat("==============================================================\n\n")
  
  cat("Economic Weights Statistics:\n")
  cat(sprintf("  Mean weight:         %.4f\n", mean(object$w)))
  cat(sprintf("  SD of weights:       %.4f\n", sd(object$w)))
  cat(sprintf("  Range:               [%.4f, %.4f]\n", min(object$w), max(object$w)))
  
  cat("\nResponse Statistics:\n")
  cat(sprintf("  Mean DeltaG:             %.4f\n", mean(object$Delta_G)))
  cat(sprintf("  SD of DeltaG:            %.4f\n", sd(object$Delta_G)))
  cat(sprintf("  Range DeltaG:            [%.4f, %.4f]\n",
              min(object$Delta_G), max(object$Delta_G)))
  
  if (!is.null(object$lpsi_comparison)) {
    cat("\nLPSI vs Base Index Comparison:\n")
    cat(sprintf("  GA improvement:      %.2f%% gain if using LPSI\n",
                (1/object$lpsi_comparison$efficiency_ratio - 1) * 100))
    cat(sprintf("  rHI improvement:     %.4f (LPSI) vs %.4f (Base)\n",
                object$lpsi_comparison$rHI_lpsi, object$rHI))
    
    # Calculate correlation between responses
    cor_responses <- cor(object$Delta_G, object$lpsi_comparison$Delta_G_lpsi)
    cat(sprintf("  Response correlation: %.4f\n", cor_responses))
    
    if (cor_responses < 0.8) {
      cat("\n[!] Low correlation between Base Index and LPSI responses.\n")
      cat("  The two methods prioritize traits differently.\n")
    }
  }
  
  cat("\n")
  invisible(object)
}
