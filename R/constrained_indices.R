#' Constrained selection indices
#' @name constrained_indices
#'
#' @description
#' Implements constrained selection index methods using the same math primitives
#' as the unrestricted Smith-Hazel index.
#'
#' Methods included:
#' - Restricted Linear Phenotypic Selection Index (RLPSI)
#' - Predetermined Proportional Gains (PPG-LPSI)
#' - Desired Gains (DG-LPSI)
#'
#' @keywords internal
NULL

.solve_sym_multi <- function(A, B) {
  B <- as.matrix(B)
  n_col <- ncol(B)
  out <- matrix(0, nrow = nrow(A), ncol = n_col)
  for (j in seq_len(n_col)) {
    out[, j] <- cpp_symmetric_solve(A, B[, j])
  }
  out
}

.index_metrics <- function(b, P, G, w = NULL, const_factor = 2.063, GAY = NULL) {
  b <- as.numeric(b)
  bPb <- cpp_quadratic_form_sym(b, P)
  bGb <- cpp_quadratic_form_sym(b, G)
  sigma_I <- if (bPb > 0) sqrt(bPb) else NA_real_
  delta_g_scalar <- if (!is.na(sigma_I)) const_factor * sigma_I else NA_real_
  delta_g_vec <- if (!is.na(sigma_I)) const_factor * (G %*% b) / sigma_I else rep(NA_real_, nrow(G))

  hI2 <- if (!is.na(bPb) && bPb > 0) bGb / bPb else NA_real_
  rHI <- if (!is.na(hI2) && hI2 >= 0) sqrt(hI2) else NA_real_

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

#' Restricted Linear Phenotypic Selection Index (RLPSI)
#'
#' @param pmat Phenotypic variance-covariance matrix
#' @param gmat Genotypic variance-covariance matrix
#' @param wmat Weight matrix
#' @param wcol Weight column number (default: 1)
#' @param restricted_traits Vector of trait indices to restrict (default: NULL).
#'   If provided, a constraint matrix C is auto-generated to enforce zero gain on these traits.
#'   Example: c(1, 3) restricts traits 1 and 3 to zero gain.
#' @param C Constraint matrix (n_traits x n_constraints). Each column is a restriction.
#'   Alternative to restricted_traits for custom constraints. Ignored if restricted_traits is provided.
#' @param GAY Genetic advance of comparative trait (optional)
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients (b.*), GA, PRE, Delta_G, rHI, hI2
#'     \item \code{b} - Numeric vector of selection index coefficients
#'     \item \code{Delta_G} - Named vector of realized correlated responses per trait
#'     \item \code{C} - Constraint matrix used
#'   }
#' @export
#' @importFrom stats setNames
#' @examples
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' wmat <- weight_mat(weight)
#' 
#' # Easy way: Restrict traits 1 and 3 to zero gain
#' result <- rlpsi(pmat, gmat, wmat, wcol = 1, restricted_traits = c(1, 3))
#' 
#' # Advanced way: Provide custom constraint matrix
#' C <- diag(ncol(pmat))[, 1, drop = FALSE]
#' result <- rlpsi(pmat, gmat, wmat, wcol = 1, C = C)
rlpsi <- function(pmat, gmat, wmat, wcol = 1, restricted_traits = NULL, C = NULL, GAY) {
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  wmat <- as.matrix(wmat)

  # Auto-create C from restricted_traits (user-friendly)
  if (!is.null(restricted_traits)) {
    if (!is.numeric(restricted_traits) || any(restricted_traits < 1) || any(restricted_traits > nrow(pmat))) {
      stop("restricted_traits must be a numeric vector of valid trait indices (1 to ", nrow(pmat), ").")
    }
    C <- diag(nrow(pmat))[, restricted_traits, drop = FALSE]
  } else if (is.null(C)) {
    stop("Either 'restricted_traits' or 'C' must be provided for RLPSI constraints.")
  }

  C <- as.matrix(C)
  if (nrow(C) != nrow(pmat)) {
    stop("C must have the same number of rows as pmat (one row per trait).")
  }

  w <- cpp_extract_vector(wmat, seq_len(nrow(pmat)), wcol - 1L)

  P_inv_G <- .solve_sym_multi(pmat, gmat)
  P_inv_Gw <- cpp_symmetric_solve(pmat, gmat %*% w)

  middle <- t(C) %*% gmat %*% P_inv_G %*% C
  middle_inv <- solve(middle)

  proj <- diag(nrow(pmat)) - P_inv_G %*% C %*% middle_inv %*% t(C) %*% gmat
  b <- proj %*% P_inv_Gw

  metrics <- .index_metrics(b, pmat, gmat, w = w, GAY = if (missing(GAY)) NULL else GAY)

  # Ensure b is a clean numeric vector (not matrix)
  b_vec <- as.numeric(b)
  b_vec <- round(b_vec, 4)
  
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))

  summary_df <- data.frame(
    b_df,
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    Delta_G = round(metrics$Delta_G, 4),
    rHI = round(metrics$rHI, 4),
    hI2 = round(metrics$hI2, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(summary_df) <- NULL

  list(
    summary = summary_df,
    b = b_vec,
    Delta_G = setNames(metrics$Delta_G_vec, colnames(pmat)),
    C = C
  )
}

#' Predetermined Proportional Gains (PPG-LPSI)
#'
#' @param pmat Phenotypic variance-covariance matrix
#' @param gmat Genotypic variance-covariance matrix
#' @param k Vector of desired proportional gains
#' @param wmat Optional weight matrix for GA/PRE calculation
#' @param wcol Weight column number (default: 1)
#' @param GAY Genetic advance of comparative trait (optional)
#'
#' @return List with summary data frame, coefficient vector, Delta_G vector, and phi
#' @export
#'
#' @examples
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' k <- rep(1, ncol(pmat))
#' ppg_lpsi(pmat, gmat, k)
ppg_lpsi <- function(pmat, gmat, k, wmat = NULL, wcol = 1, GAY) {
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  k <- as.numeric(k)

  if (length(k) != nrow(pmat)) {
    stop("k must have the same length as the number of traits.")
  }

  P_inv_G <- .solve_sym_multi(pmat, gmat)
  S <- t(gmat) %*% P_inv_G
  
  # Check for singular matrix
  S_cond <- tryCatch({
    kappa(S, exact = TRUE)
  }, error = function(e) Inf)
  
  if (is.infinite(S_cond) || S_cond > 1e10) {
    stop("Singular matrix detected in PPG-LPSI: G'P^{-1}G is rank deficient. ",
         "This can occur when traits are linearly dependent or G is not full rank.")
  }
  
  x <- solve(S, k)
  b <- P_inv_G %*% x

  w <- NULL
  if (!is.null(wmat)) {
    w <- cpp_extract_vector(as.matrix(wmat), seq_len(nrow(pmat)), wcol - 1L)
  }

  metrics <- .index_metrics(b, pmat, gmat, w = w, GAY = if (missing(GAY)) NULL else GAY)

  ratios <- metrics$Delta_G_vec / k
  phi <- if (all(is.finite(ratios) & k != 0)) mean(ratios[k != 0]) else NA_real_

  # Ensure b is a clean numeric vector (not matrix)
  b_vec <- as.numeric(b)
  b_vec <- round(b_vec, 4)
  
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))

  summary_df <- data.frame(
    b_df,
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    Delta_G = round(metrics$Delta_G, 4),
    rHI = round(metrics$rHI, 4),
    hI2 = round(metrics$hI2, 4),
    phi = round(phi, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(summary_df) <- NULL

  list(
    summary = summary_df,
    b = b_vec,
    Delta_G = setNames(metrics$Delta_G_vec, colnames(pmat)),
    phi = phi
  )
}

#' Desired Gains Index (DG-LPSI) with Implied Economic Weights
#'
#' @description
#' Implements the Pesek & Baker (1969) Desired Gains Index where breeders specify
#' target genetic gains instead of economic weights. This enhanced version includes
#' calculation of implied economic weights and feasibility checking.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param d Vector of desired genetic gains (length n_traits)
#' @param wmat (Deprecated) Not used in DG-LPSI as desired gains replace economic weights
#' @param wcol (Deprecated) Not used in DG-LPSI
#' @param GAY (Deprecated) Not used in DG-LPSI as GA/PRE are not applicable without economic weights
#' @param return_implied_weights Logical - calculate implied economic weights? (default: TRUE)
#' @param check_feasibility Logical - warn if desired gains are unrealistic? (default: TRUE)
#' @param selection_intensity Selection intensity i (default: 2.063)
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients, metrics, and implied weights
#'     \item \code{b} - Vector of selection index coefficients
#'     \item \code{Delta_G} - Named vector of achieved genetic gains per trait
#'     \item \code{desired_gains} - Named vector of desired gains (input d)
#'     \item \code{gain_errors} - Difference between desired and achieved gains
#'     \item \code{implied_weights} - Economic weights that would achieve these gains in Smith-Hazel LPSI
#'     \item \code{implied_weights_normalized} - Normalized implied weights (max absolute = 1)
#'     \item \code{feasibility} - Data frame with feasibility analysis per trait
#'     \item \code{hI2} - Index heritability
#'     \item \code{rHI} - Index accuracy
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#' 
#' 1. Index coefficients: \eqn{\mathbf{b} = \mathbf{G}^{-1} \mathbf{d}}
#' 
#' 2. Expected response: \eqn{\Delta \mathbf{G} = \mathbf{G}\mathbf{b}}
#' 
#' 3. Implied economic weights (Section 1.4 of Chapter 4): 
#'    \deqn{\hat{\mathbf{w}} = \mathbf{G}^{-1} \mathbf{P} \mathbf{b}}
#' 
#' The implied weights represent the economic values that would have been needed
#' in a Smith-Hazel index to achieve the desired gains. Large implied weights
#' indicate traits that are "expensive" to improve (low heritability or unfavorable
#' correlations), while small weights indicate traits that are "cheap" to improve.
#'
#' \strong{Feasibility Checking:}
#' 
#' The function estimates maximum possible gains as approximately 3.0 * sqrt(G_ii)
#' (assuming very intense selection with i ≈ 3.0) and warns if desired gains
#' exceed 80% of these theoretical maxima.
#'
#' @references
#' Pesek, J., & Baker, R. J. (1969). Desired improvement in relation to 
#' selection indices. \emph{Canadian Journal of Plant Science}, 49(6), 803-804.
#'
#' @export
#'
#' @examples
#' # Load data
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' 
#' # Specify desired gains (e.g., increase each trait by 1 unit)
#' desired_gains <- rep(1, ncol(gmat))
#' 
#' # Calculate Desired Gains Index with all enhancements
#' result <- dg_lpsi(pmat, gmat, d = desired_gains)
#' 
#' # View summary
#' print(result$summary)
#' 
#' # Extract implied weights to understand relative "cost" of gains
#' print(result$implied_weights_normalized)
#' 
#' # Check feasibility
#' print(result$feasibility)
dg_lpsi <- function(pmat, gmat, d, wmat = NULL, wcol = 1, GAY,
                    return_implied_weights = TRUE,
                    check_feasibility = TRUE,
                    selection_intensity = 2.063) {
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================
  
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  d <- as.numeric(d)
  n_traits <- nrow(pmat)

  if (length(d) != n_traits) {
    stop("Length of d must equal number of traits (nrow(pmat) = ", n_traits, ").")
  }
  
  if (nrow(pmat) != ncol(pmat) || nrow(gmat) != ncol(gmat)) {
    stop("pmat and gmat must be square matrices.")
  }
  
  if (nrow(pmat) != nrow(gmat)) {
    stop("pmat and gmat must have the same dimensions.")
  }

  # ============================================================================
  # STEP 1: Calculate Index Coefficients (Section 1.2)
  # Formula: b = G^(-1) d
  # ============================================================================

  P_inv_G <- .solve_sym_multi(pmat, gmat)
  S <- t(gmat) %*% P_inv_G
  
  # Check for singular matrix
  S_cond <- tryCatch({
    kappa(S, exact = TRUE)
  }, error = function(e) Inf)
  
  if (is.infinite(S_cond) || S_cond > 1e10) {
    stop("Singular matrix detected in DG-LPSI: G'P^{-1}G is rank deficient. ",
         "This can occur when traits are linearly dependent or G is not full rank. ",
         "Consider using a subset of traits or checking your variance-covariance matrices.")
  }
  
  x <- solve(S, d)
  b <- P_inv_G %*% x
  b <- as.numeric(b)
  
  # Check for numerical issues
  if (any(is.na(b)) || any(is.infinite(b))) {
    stop("Index coefficients contain NA or Inf. Check that gmat is invertible.")
  }

  # ============================================================================
  # STEP 2: Calculate Expected Response (Section 1.3)
  # Formula: ΔG = G * b
  # ============================================================================
  
  Delta_G_vec <- gmat %*% b
  Delta_G_vec <- as.numeric(Delta_G_vec)
  
  # Calculate error between desired and achieved gains
  gain_errors <- Delta_G_vec - d
  max_error <- max(abs(gain_errors))
  
  if (max_error > 1e-4) {
    warning(
      "Desired gains not perfectly achieved. Maximum error: ", 
      format(max_error, scientific = TRUE, digits = 4),
      "\nThis may indicate numerical instability or rank-deficient gmat."
    )
  }

  # ============================================================================
  # STEP 3: Calculate Standard Metrics
  # ============================================================================
  
  metrics <- .index_metrics(b, pmat, gmat, w = NULL, 
                            const_factor = selection_intensity, 
                            GAY = NULL)

  # ============================================================================
  # STEP 4: Calculate Implied Economic Weights (Section 1.4) ⭐ NEW
  # Formula: ŵ = G^(-1) P b
  # 
  # Interpretation: These are the economic weights that would have been 
  # needed in a Smith-Hazel index to achieve the desired gains.
  # ============================================================================
  
  implied_weights <- NULL
  implied_weights_normalized <- NULL
  
  if (return_implied_weights) {
    # Use generalized inverse for numerical stability
    gmat_inv <- MASS::ginv(gmat)
    implied_weights <- gmat_inv %*% pmat %*% b
    implied_weights <- as.numeric(implied_weights)
    
    # Check for numerical issues
    if (any(is.na(implied_weights)) || any(is.infinite(implied_weights))) {
      warning("Implied weights contain NA or Inf. Check matrix conditioning.")
      implied_weights <- rep(NA_real_, n_traits)
      implied_weights_normalized <- rep(NA_real_, n_traits)
    } else {
      # Normalize for interpretability (largest absolute weight = 1)
      max_abs_weight <- max(abs(implied_weights))
      
      if (max_abs_weight > 0) {
        implied_weights_normalized <- implied_weights / max_abs_weight
      } else {
        implied_weights_normalized <- implied_weights
      }
      
      # Name vectors
      trait_names <- colnames(pmat)
      if (!is.null(trait_names)) {
        names(implied_weights) <- trait_names
        names(implied_weights_normalized) <- trait_names
      }
    }
  }

  # ============================================================================
  # STEP 5: Feasibility Check ⭐ NEW
  # Check if desired gains are realistic given genetic variances
  # ============================================================================
  
  feasibility_metrics <- NULL
  
  if (check_feasibility) {
    # Maximum possible gains under very intense selection
    # Approximation: ΔG_max ≈ i * sqrt(G_ii) using provided selection intensity
    genetic_sd <- sqrt(diag(gmat))
    max_possible_gains <- selection_intensity * genetic_sd
    
    # Ratio of desired to maximum possible
    gain_ratios <- abs(d) / max_possible_gains
    
    # Flag unrealistic gains (> 80% of theoretical maximum)
    unrealistic_traits <- which(gain_ratios > 0.8)
    
    if (length(unrealistic_traits) > 0) {
      trait_names_warn <- if (!is.null(colnames(gmat))) {
        colnames(gmat)[unrealistic_traits]
      } else {
        paste0("Trait_", unrealistic_traits)
      }
      
      warning(
        "Desired gains may be unrealistic for trait(s): ",
        paste(trait_names_warn, collapse = ", "),
        "\nDesired gains exceed 80% of theoretical maximum (i = ", 
        round(selection_intensity, 2), ").",
        "\nConsider reducing targets or using higher selection intensity.",
        call. = FALSE
      )
    }
    
    # Build feasibility data frame
    trait_names <- colnames(gmat)
    if (is.null(trait_names)) {
      trait_names <- paste0("Trait_", seq_len(n_traits))
    }
    
    feasibility_metrics <- data.frame(
      trait = trait_names,
      desired_gain = d,
      achieved_gain = Delta_G_vec,
      genetic_sd = genetic_sd,
      max_possible_gain = round(max_possible_gains, 4),
      feasibility_ratio = round(gain_ratios, 4),
      is_realistic = gain_ratios <= 0.8,
      stringsAsFactors = FALSE
    )
  }

  # ============================================================================
  # STEP 6: Build Summary Data Frame
  # ============================================================================
  
  b_vec <- round(b, 4)
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))
  
  summary_df <- data.frame(
    b_df,
    Delta_G = round(metrics$Delta_G, 4),
    hI2 = round(metrics$hI2, 4),
    rHI = round(metrics$rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Add implied weights to summary if calculated
  if (return_implied_weights && !is.null(implied_weights)) {
    for (i in seq_along(implied_weights)) {
      summary_df[[paste0("implied_w.", i)]] <- round(implied_weights[i], 4)
    }
    for (i in seq_along(implied_weights_normalized)) {
      summary_df[[paste0("implied_w_norm.", i)]] <- round(implied_weights_normalized[i], 4)
    }
  }
  
  rownames(summary_df) <- NULL

  # ============================================================================
  # STEP 7: Return Enhanced Results
  # ============================================================================
  
  result <- list(
    summary = summary_df,
    b = b_vec,
    Delta_G = setNames(Delta_G_vec, colnames(pmat)),
    desired_gains = setNames(d, colnames(pmat)),
    gain_errors = setNames(gain_errors, colnames(pmat)),
    implied_weights = implied_weights,
    implied_weights_normalized = implied_weights_normalized,
    feasibility = feasibility_metrics,
    hI2 = metrics$hI2,
    rHI = metrics$rHI,
    selection_intensity = selection_intensity,
    method = "Desired Gains Index (Pesek & Baker, 1969)"
  )
  
  class(result) <- c("dg_lpsi", "selection_index", "list")
  return(result)
}

#' Print method for Desired Gains Index
#'
#' @param x Object of class 'dg_lpsi'
#' @param ... Additional arguments (unused)
#' @export
print.dg_lpsi <- function(x, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════════════════\n")
  cat("DESIRED GAINS INDEX (Pesek & Baker, 1969)\n")
  cat("══════════════════════════════════════════════════════════════\n\n")
  
  cat("Selection intensity (i):", x$selection_intensity, "\n\n")
  
  cat("Index Coefficients (b = G^-1 d):\n")
  print(round(x$b, 4))
  
  cat("\n\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat("DESIRED vs. ACHIEVED GAINS\n")
  cat("─────────────────────────────────────────────────────────────\n")
  
  # Handle missing trait names
  trait_names <- names(x$desired_gains)
  if (is.null(trait_names) || length(trait_names) == 0) {
    trait_names <- paste0("Trait_", seq_along(x$desired_gains))
  }
  
  comparison <- data.frame(
    Trait = trait_names,
    Desired = round(as.numeric(x$desired_gains), 4),
    Achieved = round(as.numeric(x$Delta_G), 4),
    Error = round(as.numeric(x$gain_errors), 4),
    stringsAsFactors = FALSE
  )
  print(comparison)
  
  max_error <- max(abs(x$gain_errors))
  if (max_error < 1e-4) {
    cat("\n✓ Desired gains achieved with high precision\n")
  } else if (max_error < 0.01) {
    cat("\n⚠ Small errors present (max:", format(max_error, scientific = TRUE, digits = 4), ")\n")
  } else {
    cat("\n✗ Significant errors detected (max:", round(max_error, 6), ")\n")
    cat("  Check for numerical instability or rank-deficient gmat\n")
  }
  
  cat("\n\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat("INDEX METRICS\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat("├─ Index Heritability (h²_I):", round(x$hI2, 4), "\n")
  cat("└─ Index Accuracy (r_HI):", round(x$rHI, 4), "\n")
  
  if (!is.null(x$implied_weights)) {
    cat("\n\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat("IMPLIED ECONOMIC WEIGHTS (ŵ = G^-1 P b)\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat("These are the economic weights that would have been needed in\n")
    cat("a Smith-Hazel index to achieve the desired gains.\n\n")
    
    # Handle missing trait names
    weight_names <- names(x$implied_weights)
    if (is.null(weight_names) || length(weight_names) == 0) {
      weight_names <- paste0("Trait_", seq_along(x$implied_weights))
    }
    
    weights_df <- data.frame(
      Trait = weight_names,
      Implied_Weight = round(as.numeric(x$implied_weights), 4),
      Normalized = round(as.numeric(x$implied_weights_normalized), 4),
      stringsAsFactors = FALSE
    )
    print(weights_df)
    
    cat("\n📊 Interpretation:\n")
    cat("   • Large weights → trait is 'expensive' to improve\n")
    cat("   • Small weights → trait is 'cheap' to improve\n")
  }
  
  if (!is.null(x$feasibility)) {
    cat("\n\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat("FEASIBILITY ANALYSIS\n")
    cat("─────────────────────────────────────────────────────────────\n")
    print(x$feasibility)
    
    n_unrealistic <- sum(!x$feasibility$is_realistic)
    if (n_unrealistic > 0) {
      cat("\n⚠ WARNING:", n_unrealistic, "trait(s) have unrealistic desired gains\n")
      cat("  Gains exceed 80% of theoretical maximum (i = 3.0)\n")
      cat("  Consider reducing targets or increasing selection intensity.\n")
    } else {
      cat("\n✓ All desired gains are feasible\n")
    }
  }
  
  cat("\n")
  cat("══════════════════════════════════════════════════════════════\n")
  cat("\n")
  invisible(x)
}

#' Summary method for Desired Gains Index
#'
#' @param object Object of class 'dg_lpsi'
#' @param ... Additional arguments (unused)
#' @export
summary.dg_lpsi <- function(object, ...) {
  print(object, ...)
  
  cat("\n📚 ADDITIONAL DETAILS:\n\n")
  
  cat("Summary Statistics:\n")
  cat("├─ Mean desired gain:", round(mean(object$desired_gains), 4), "\n")
  cat("├─ Mean achieved gain:", round(mean(object$Delta_G), 4), "\n")
  cat("├─ Mean absolute error:", format(mean(abs(object$gain_errors)), scientific = TRUE, digits = 4), "\n")
  
  if (!is.null(object$implied_weights_normalized)) {
    cat("└─ Mean implied weight (normalized):", 
        round(mean(object$implied_weights_normalized), 4), "\n")
  }
  
  cat("\n")
  invisible(object)
}
