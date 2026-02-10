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

#' Desired Gains Index (DG-LPSI)
#'
#' @param pmat Phenotypic variance-covariance matrix
#' @param gmat Genotypic variance-covariance matrix
#' @param d Vector of desired gains
#' @param wmat (Deprecated) Not used in DG-LPSI as desired gains replace economic weights
#' @param wcol (Deprecated) Not used in DG-LPSI
#' @param GAY (Deprecated) Not used in DG-LPSI as GA/PRE are not applicable without economic weights
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients (b.*), Delta_G, hI2, and rHI
#'     \item \code{b} - Vector of selection index coefficients
#'     \item \code{Delta_G} - Vector of realized correlated responses per trait
#'     \item \code{desired_gains} - Vector of desired gains (input d)
#'   }
#'
#' @details
#' The Desired Gains Index replaces economic weights with desired gains.
#' GA (Genetic Advance) and PRE (Percent Relative Efficiency) are not calculated
#' as they require economic weights. Only index heritability (hI2) and accuracy (rHI)
#' are meaningful metrics for DG-LPSI.
#'
#' @export
#'
#' @examples
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' d <- rep(1, ncol(pmat))
#' result <- dg_lpsi(pmat, gmat, d)
#' print(result$summary)
dg_lpsi <- function(pmat, gmat, d, wmat = NULL, wcol = 1, GAY) {
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  d <- as.numeric(d)

  if (length(d) != nrow(pmat)) {
    stop("d must have the same length as the number of traits.")
  }

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

  # DG-LPSI replaces economic weights with desired gains
  # GA and PRE are not applicable without economic weights
  # Only calculate index heritability (hI2) and related metrics
  metrics <- .index_metrics(b, pmat, gmat, w = NULL, GAY = NULL)

  # Ensure b is a clean numeric vector (not matrix)
  b_vec <- as.numeric(b)
  b_vec <- round(b_vec, 4)
  
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))

  # Note: GA and PRE are NA for DG-LPSI (no economic weights)
  # Only hI2 (index heritability) and rHI (accuracy of index) are meaningful
  summary_df <- data.frame(
    b_df,
    Delta_G = round(metrics$Delta_G, 4),
    hI2 = round(metrics$hI2, 4),
    rHI = round(metrics$rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(summary_df) <- NULL

  list(
    summary = summary_df,
    b = b_vec,
    Delta_G = setNames(metrics$Delta_G_vec, colnames(pmat)),
    desired_gains = setNames(d, colnames(pmat))
  )
}
