#' Constrained Phenotypic Selection Indices (Chapter 3)
#' @name constrained_indices
#'
#' @description
#' Implements constrained phenotypic selection index methods from Chapter 3.
#' These methods allow breeders to impose restrictions on genetic gains
#' or specify target gains while maintaining selection efficiency.
#'
#' Methods included:
#' - Restricted Linear Phenotypic Selection Index (RLPSI) - Kempthorne & Nordskog (1959)
#' - Predetermined Proportional Gains (PPG-LPSI) - Tallis (1962)
#' - Desired Gains Index (DG-LPSI) - Pesek & Baker (1969)
#'
#' All implementations use C++ primitives for mathematical operations.
#'
#' @references
#' Kempthorne, O., & Nordskog, A. W. (1959). Restricted selection indices.
#' Biometrics, 15(1), 10-19.
#'
#' Tallis, G. M. (1962). A selection index for optimum genotype.
#' Biometrics, 18(1), 120-122.
#'
#' Pesek, J., & Baker, R. J. (1969). Desired improvement in relation to selection indices.
#' Canadian Journal of Plant Science, 49(6), 803-804.
#'
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 3.
#'
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom MASS ginv
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
#' @description
#' Implements the Restricted LPSI where genetic gains are constrained to zero
#' for specific traits while maximizing gains for others.
#' Based on Kempthorne & Nordskog (1959).
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param wmat Weight matrix (n_traits x k), or vector
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
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 3, Section 3.1):}
#'
#' The RLPSI minimizes the mean squared difference between I = b'y and H = w'g
#' subject to the restriction: C'ΔG = 0
#'
#' Coefficient formula:
#' \deqn{b_r = [I - P^{-1}GC(C'GP^{-1}GC)^{-1}C'G]P^{-1}Gw}
#'
#' Where:
#' - P = Phenotypic variance-covariance matrix
#' - G = Genotypic variance-covariance matrix
#' - C = Constraint matrix (each column enforces one restriction)
#' - w = Economic weights
#'
#' The constraint C'ΔG = 0 ensures zero genetic gain for restricted traits.
#'
#' @export
#' @examples
#' \dontrun{
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
#' }
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

  # Use ginv for robustness when restricting collinear traits
  middle <- t(C) %*% gmat %*% P_inv_G %*% C
  middle_inv <- ginv(middle)  # Robust to singular matrices

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
#' @description
#' Implements the PPG-LPSI where breeders specify desired proportional gains
#' between traits rather than restricting specific traits to zero.
#' Based on Tallis (1962).
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param k Vector of desired proportional gains (length n_traits).
#'   Example: k = c(2, 1, 1) means trait 1 should gain twice as much as traits 2 and 3.
#' @param wmat Optional weight matrix for GA/PRE calculation
#' @param wcol Weight column number (default: 1)
#' @param GAY Genetic advance of comparative trait (optional)
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients and metrics
#'     \item \code{b} - Vector of PPG-LPSI coefficients
#'     \item \code{Delta_G} - Expected genetic gains per trait
#'     \item \code{phi} - Proportionality constant
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 3, Section 3.2):}
#'
#' The PPG-LPSI achieves gains in specific proportions: ΔG = φk
#'
#' Coefficient formula (Tallis, 1962):
#' \\deqn{b = P^{-1}G(G'P^{-1}G)^{-1}k}
#'
#' Where:
#' - k = Vector of desired proportions
#' - φ = Proportionality constant (determined by selection intensity and variances)
#'
#' The constraint ensures ΔG₁:ΔG₂:ΔG₃ = k₁:k₂:k₃
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' 
#' # Gains in ratio 2:1:1:1:1:1:1
#' k <- c(2, 1, 1, 1, 1, 1, 1)
#' result <- ppg_lpsi(pmat, gmat, k)
#' }
ppg_lpsi <- function(pmat, gmat, k, wmat = NULL, wcol = 1, GAY) {
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  k <- as.numeric(k)

  if (length(k) != nrow(pmat)) {
    stop("k must have the same length as the number of traits.")
  }

  # Correct PPG-LPSI formula (Tallis, 1962): b = P^{-1}GP^{-1}k
  # WARNING: Do NOT use P^{-1}G(G'P^{-1}G)^{-1}k as it simplifies to G^{-1}k
  P_inv_G <- .solve_sym_multi(pmat, gmat)  # P^{-1}G
  P_inv_k <- cpp_symmetric_solve(pmat, k)   # P^{-1}k
  b <- P_inv_G %*% P_inv_k                   # P^{-1}G(P^{-1}k) = P^{-1}GP^{-1}k
  
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

#' Desired Gains Index (DG-LPSI)
#'
#' @description
#' Implements the Pesek & Baker (1969) Desired Gains Index where breeders specify
#' target genetic gains instead of economic weights. This enhanced version includes
#' calculation of implied economic weights and feasibility checking.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param d Vector of desired genetic gains (length n_traits).
#'   Example: d = c(1.5, 0.8, -0.2) means gain +1.5 in trait 1, +0.8 in trait 2, -0.2 in trait 3.
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
#' 1. Index coefficients: \eqn{\mathbf{b} = \mathbf{G}^{-1}\mathbf{d}}
#'
#' 2. Expected response: \eqn{\Delta \mathbf{G} = (i/\sigma_I) \mathbf{G}\mathbf{b}}
#'
#' \strong{CRITICAL: Scale Invariance Property}
#'
#' The achieved gains \eqn{\Delta\mathbf{G}} are determined by selection intensity (i),
#' genetic variance (G), and phenotypic variance (P), NOT by scaling \eqn{\mathbf{b}}.
#' If you multiply \eqn{\mathbf{b}} by constant c, \eqn{\sigma_I} also scales by c, causing
#' complete cancellation in \eqn{\Delta\mathbf{G} = (i/(c\sigma_I))\mathbf{G}(c\mathbf{b}) = (i/\sigma_I)\mathbf{G}\mathbf{b}}.
#'
#' \strong{What DG-LPSI Actually Achieves:}
#'
#' - Proportional gains matching the RATIOS in d (not absolute magnitudes)
#' - Achieved magnitude depends on biological/genetic constraints
#' - Use feasibility checking to verify if desired gains are realistic
#'
#' 3. Implied economic weights (Section 1.4 of Chapter 4):
#'    \deqn{\hat{\mathbf{w}} = \mathbf{G}^{-1} \mathbf{P} \mathbf{b}}
#'
#' The implied weights represent the economic values that would have been needed
#' in a Smith-Hazel index to achieve the desired gain PROPORTIONS. Large implied weights
#' indicate traits that are "expensive" to improve (low heritability or unfavorable
#' correlations), while small weights indicate traits that are "cheap" to improve.
#'
#' \strong{Feasibility Checking:}
#'
#' The function estimates maximum possible gains as approximately 3.0 * sqrt(G_ii)
#' (assuming very intense selection with i ~ 3.0) and warns if desired gains
#' exceed 80% of these theoretical maxima.
#'
#' @references
#' Pesek, J., & Baker, R. J. (1969). Desired improvement in relation to
#' selection indices. \emph{Canadian Journal of Plant Science}, 49(6), 803-804.
#'
#' @importFrom MASS ginv
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
  # Formula: b = G^(-1) d (Pesek & Baker, 1969)
  # ============================================================================

  # Use ginv for robustness (handles near-singular cases)
  gmat_inv <- ginv(gmat)
  b <- gmat_inv %*% d
  b <- as.numeric(b)

  # Check for numerical issues
  if (any(is.na(b)) || any(is.infinite(b))) {
    stop("Index coefficients contain NA or Inf. Check that gmat is invertible.")
  }

  # ============================================================================
  # STEP 2: Calculate Expected Response (Section 1.3)
  # Formula: DeltaG = (i/sigma_I) * G * b
  # CRITICAL: Must use .index_metrics to get correctly scaled Delta_G
  # ============================================================================

  # ============================================================================
  # STEP 3: Calculate Standard Metrics (includes correct Delta_G calculation)
  # ============================================================================

  metrics <- .index_metrics(b, pmat, gmat, w = NULL,
                            const_factor = selection_intensity,
                            GAY = NULL)
  
  # Extract correctly scaled achieved gains from metrics
  Delta_G_vec <- metrics$Delta_G_vec
  
  # Calculate proportional match (Pesek & Baker only guarantees proportions)
  # Check if achieved gains are proportional to desired gains
  gain_ratios <- Delta_G_vec / d
  gain_ratios[abs(d) < 1e-10] <- NA_real_  # Avoid division by zero
  
  # Check if ratios are consistent (all approximately equal)
  valid_ratios <- gain_ratios[is.finite(gain_ratios)]
  if (length(valid_ratios) > 1) {
    ratio_consistency <- sd(valid_ratios) / mean(valid_ratios)
    if (ratio_consistency > 0.01) {
      warning(
        "Achieved gains are not perfectly proportional to desired gains.\n",
        "Coefficient of variation in ratios: ", round(ratio_consistency * 100, 2), "%\n",
        "This may indicate numerical instability or ill-conditioned matrices.",
        call. = FALSE
      )
    }
  }
  
  # Calculate proportional scale factor
  avg_ratio <- if (length(valid_ratios) > 0) mean(valid_ratios) else 1
  
  # "Error" measures deviation from perfect proportionality (not absolute difference)
  gain_errors <- Delta_G_vec - (avg_ratio * d)

  # ============================================================================
  # STEP 4: Calculate Implied Economic Weights (Section 1.4) [NEW]
  # Formula: w-hat = G^(-1) P b
  #
  # Interpretation: These are the economic weights that would have been
  # needed in a Smith-Hazel index to achieve the desired gain PROPORTIONS.
  # ============================================================================

  implied_weights <- NULL
  implied_weights_normalized <- NULL

  if (return_implied_weights) {
    # gmat_inv already computed in STEP 1
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
  # STEP 5: Feasibility Check [NEW]
  # Check if desired gains are realistic given genetic variances
  # ============================================================================

  feasibility_metrics <- NULL

  if (check_feasibility) {
    # Maximum possible gains under very intense selection
    # Approximation: DeltaG_max ~ i * sqrt(G_ii) using provided selection intensity
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
  cat("==============================================================\n")
  cat("DESIRED GAINS INDEX (Pesek & Baker, 1969)\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (i):", x$selection_intensity, "\n\n")

  cat("Index Coefficients (b = G^-1 d):\n")
  print(round(x$b, 4))

  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("DESIRED vs. ACHIEVED GAINS\n")
  cat("-------------------------------------------------------------\n")

  # Handle missing trait names
  trait_names <- names(x$desired_gains)
  if (is.null(trait_names) || length(trait_names) == 0) {
    trait_names <- paste0("Trait_", seq_along(x$desired_gains))
  }

  # Calculate proportional scale and ratios
  gain_ratios <- as.numeric(x$Delta_G) / as.numeric(x$desired_gains)
  gain_ratios[abs(as.numeric(x$desired_gains)) < 1e-10] <- NA
  avg_scale <- mean(gain_ratios, na.rm = TRUE)
  
  comparison <- data.frame(
    Trait = trait_names,
    Desired = round(as.numeric(x$desired_gains), 4),
    Achieved = round(as.numeric(x$Delta_G), 4),
    Proportional = round(gain_ratios, 4),
    Prop_Error = round(as.numeric(x$gain_errors), 6),
    stringsAsFactors = FALSE
  )
  print(comparison)

  cat("\n** CRITICAL NOTE: Pesek & Baker guarantees PROPORTIONAL gains, not absolute magnitudes.\n")
  cat("   Average proportional scale factor (φ):", round(avg_scale, 4), "\n")
  cat("   If all 'Proportional' values are equal, proportionality is achieved.\n")
  
  max_error <- max(abs(x$gain_errors))
  if (max_error < 1e-4) {
    cat("\n[OK] Desired gain proportions achieved with high precision\n")
  } else if (max_error < 0.01) {
    cat("\n[WARNING] Small proportionality errors present (max:", format(max_error, scientific = TRUE, digits = 4), ")\n")
  } else {
    cat("\n[ERROR] Significant proportionality errors detected (max:", round(max_error, 6), ")\n")
    cat("  Check for numerical instability or rank-deficient gmat\n")
  }

  cat("\n\n")
  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat("|- Index Heritability (h2_I):", round(x$hI2, 4), "\n")
  cat("|- Index Accuracy (r_HI):", round(x$rHI, 4), "\n")

  if (!is.null(x$implied_weights)) {
    cat("\n\n")
    cat("-------------------------------------------------------------\n")
    cat("IMPLIED ECONOMIC WEIGHTS (w-hat = G^-1 P b)\n")
    cat("-------------------------------------------------------------\n")
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

    cat("\n** Interpretation:\n")
    cat("   - Large weights = trait is 'expensive' to improve\n")
    cat("   - Small weights = trait is 'cheap' to improve\n")
  }

  if (!is.null(x$feasibility)) {
    cat("\n\n")
    cat("-------------------------------------------------------------\n")
    cat("FEASIBILITY ANALYSIS\n")
    cat("-------------------------------------------------------------\n")
    print(x$feasibility)

    n_unrealistic <- sum(!x$feasibility$is_realistic)
    if (n_unrealistic > 0) {
      cat("\n[!] WARNING:", n_unrealistic, "trait(s) have unrealistic desired gains\n")
      cat("  Gains exceed 80% of theoretical maximum (i = 3.0)\n")
      cat("  Consider reducing targets or increasing selection intensity.\n")
    } else {
      cat("\n[OK] All desired gains are feasible\n")
    }
  }

  cat("\n")
  cat("==============================================================\n")
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

  cat("\n[INFO] ADDITIONAL DETAILS:\n\n")

  cat("Summary Statistics:\n")
  cat("|- Mean desired gain:", round(mean(object$desired_gains), 4), "\n")
  cat("|- Mean achieved gain:", round(mean(object$Delta_G), 4), "\n")
  cat("|- Mean absolute error:", format(mean(abs(object$gain_errors)), scientific = TRUE, digits = 4), "\n")

  if (!is.null(object$implied_weights_normalized)) {
    cat("|- Mean implied weight (normalized):",
        round(mean(object$implied_weights_normalized), 4), "\n")
  }

  cat("\n")
  invisible(object)
}
