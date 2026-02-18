#' Linear Phenotypic Eigen Selection Index Methods (Chapter 7)
#' @name eigen_indices
#'
#' @description
#' Implements the Linear Phenotypic Eigen Selection Index methods from Chapter 7.
#' These methods resolve index coefficients by maximizing the accuracy squared
#' (rho_HI^2) through an eigenvalue problem rather than requiring pre-specified
#' economic weights.
#'
#' Methods included:
#' - ESIM   : Linear Phenotypic Eigen Selection Index (Section 7.1)
#' - RESIM  : Linear Phenotypic Restricted Eigen Selection Index (Section 7.2)
#' - PPG-ESIM: Predetermined Proportional Gain Eigen Selection Index (Section 7.3)
#'
#' All implementations use C++ primitives (math_primitives.cpp) for quadratic forms
#' and symmetric solves, while eigendecompositions use R's eigen() for correctness
#' and compatibility with the existing package architecture.
#'
#' @section Mathematical Foundation:
#'
#' Unlike classical LPSI which requires economic weights w, the ESIM family resolves
#' the index vector b_E by maximizing the squared accuracy:
#' \deqn{\rho_{HI}^2 = \frac{b'Cb}{b'Pb}}
#' leading to the generalized eigenproblem \eqn{(\mathbf{P}^{-1}\mathbf{C} - \lambda^2 I)b = 0}.
#'
#' The largest eigenvalue lambda_E^2 equals the maximum achievable index heritability,
#' and the corresponding eigenvector b_E contains the optimal index coefficients.
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant
#' Breeding. Springer International Publishing. Chapter 7.
#'
#' Ceron-Rojas, J. J., Crossa, J., Sahagun-Castellanos, J., Castillo-Gonzalez, F.,
#' & Santacruz-Varela, A. (2006). A selection index method based on eigen analysis.
#' Crop Science, 46(4), 1711-1721.
#'
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom MASS ginv
NULL

# ==============================================================================
# LOCAL HELPERS (mirror constrained_indices.R pattern for self-contained module)
# ==============================================================================

#' Solve symmetric linear system for multiple right-hand sides
#' @keywords internal
.esim_solve_sym_multi <- function(A, B) {
  B <- as.matrix(B)
  out <- matrix(0, nrow = nrow(A), ncol = ncol(B))
  for (j in seq_len(ncol(B))) {
    out[, j] <- cpp_symmetric_solve(A, B[, j])
  }
  out
}

#' Compute standard metrics for an eigen-based index
#'
#' @description
#' Computes all standard metrics using C++ math primitives.
#' For ESIM family indices the eigenvalue lambda^2 IS the index heritability.
#'
#' @param b  Index coefficient vector (eigenvector or transformed eigenvector)
#' @param P  Phenotypic VCV matrix
#' @param G  Genotypic VCV matrix (C in LaTeX notation)
#' @param lambda2  Eigenvalue associated with b (used for h^2_I when exact)
#' @param k_I  Selection intensity constant (default 2.063)
#' @keywords internal
.eigen_index_metrics <- function(b, P, G, lambda2 = NULL, k_I = 2.063) {
  b <- as.numeric(b)

  # --- Quadratic forms via C++ primitives ---
  bPb  <- cpp_quadratic_form_sym(b, P)   # b'Pb
  bGb  <- cpp_quadratic_form_sym(b, G)   # b'Gb

  sigma_I <- if (bPb > 0) sqrt(bPb) else NA_real_

  # Selection response: R = k_I * sigma_I
  Delta_G_scalar <- if (!is.na(sigma_I)) k_I * sigma_I else NA_real_

  # Expected genetic gain per trait: E = (k_I / sigma_I) * Gb
  Delta_G_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    as.vector(k_I * (G %*% b) / sigma_I)
  } else {
    rep(NA_real_, nrow(G))
  }

  # Index heritability: use eigenvalue if available (exact), else b'Gb / b'Pb
  hI2 <- if (!is.null(lambda2)) {
    as.numeric(lambda2)
  } else if (!is.na(bPb) && bPb > 0) {
    bGb / bPb
  } else {
    NA_real_
  }

  # Accuracy: r_HI = sqrt(h^2_I)
  rHI <- if (!is.na(hI2) && hI2 >= 0) sqrt(hI2) else NA_real_

  list(
    bPb        = bPb,
    bGb        = bGb,
    sigma_I    = sigma_I,
    Delta_G    = Delta_G_scalar,
    Delta_G_vec = Delta_G_vec,
    hI2        = hI2,
    rHI        = rHI
  )
}

#' Select the leading eigenvector from a real-eigendecomposition
#'
#' @description
#' Returns the eigenvector paired with the *largest real* eigenvalue that is
#' positive (negative eigenvalues indicate numerical noise or rank deficiency).
#' Normalises the sign so that the first non-zero element is positive.
#'
#' @param mat Square matrix to decompose
#' @param tol Eigenvalue tolerance (default 1e-8)
#' @return List: \code{vector}, \code{value}, \code{all_values}
#' @keywords internal
.leading_eigenvector <- function(mat, tol = 1e-8) {
  ev   <- eigen(mat, symmetric = FALSE)
  vals <- Re(ev$values)           # work with real parts
  vecs <- Re(ev$vectors)

  # Keep only positive eigenvalues
  pos  <- which(vals > tol)
  if (length(pos) == 0) {
    stop("No positive eigenvalues found. ",
         "Check that pmat and gmat are valid variance-covariance matrices.")
  }

  # Leading (largest positive) eigenvalue
  idx  <- pos[which.max(vals[pos])]
  bvec <- vecs[, idx]

  # Canonical sign: make the largest-magnitude element positive
  bvec <- bvec * sign(bvec[which.max(abs(bvec))])

  list(vector = bvec, value = vals[idx], all_values = vals)
}

# ==============================================================================
# 7.1  ESIM - Linear Phenotypic Eigen Selection Index
# ==============================================================================

#' Linear Phenotypic Eigen Selection Index (ESIM)
#'
#' @description
#' Implements the ESIM by maximising the squared accuracy \eqn{\rho_{HI}^2}
#' through the generalised eigenproblem of the multi-trait heritability matrix
#' \eqn{\mathbf{P}^{-1}\mathbf{C}}.
#'
#' Unlike the Smith-Hazel LPSI, **no economic weights are required**.  The net
#' genetic merit vector \eqn{\mathbf{w}_E} is instead implied by the solution.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#'   Corresponds to \strong{C} in the Chapter 7 notation.
#' @param selection_intensity Selection intensity constant \eqn{k_I}
#'   (default: 2.063 for 10\% selection).
#' @param n_indices Number of leading ESIM vectors to return (default: 1).
#'   Returning >1 provides a ranked set of indices for comparative analysis.
#'
#' @return Object of class \code{"esim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with b coefficients, hI2, rHI, sigma_I,
#'     Delta_G, and lambda2 for each index requested.}
#'   \item{\code{b}}{Named numeric vector of optimal ESIM coefficients (1st index).}
#'   \item{\code{Delta_G}}{Named numeric vector of expected genetic gains per trait.}
#'   \item{\code{sigma_I}}{Standard deviation of the index \eqn{\sigma_I}.}
#'   \item{\code{hI2}}{Index heritability \eqn{h^2_{I_E}} (= leading eigenvalue).}
#'   \item{\code{rHI}}{Accuracy \eqn{r_{HI_E}}.}
#'   \item{\code{lambda2}}{Leading eigenvalue (maximised index heritability).}
#'   \item{\code{implied_w}}{Implied economic weights \eqn{\mathbf{w}_E}.}
#'   \item{\code{all_eigenvalues}}{All eigenvalues of \eqn{\mathbf{P}^{-1}\mathbf{C}}.}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Eigenproblem (Section 7.1):}
#' \deqn{(\mathbf{P}^{-1}\mathbf{C} - \lambda_E^2 \mathbf{I})\mathbf{b}_E = 0}
#'
#' The solution \eqn{\lambda_E^2} (largest eigenvalue) equals the maximum
#' achievable index heritability \eqn{h^2_{I_E}}.
#'
#' \strong{Key metrics:}
#' \deqn{R_E = k_I \sqrt{\mathbf{b}_E^{\prime}\mathbf{P}\mathbf{b}_E}}
#' \deqn{\mathbf{E}_E = k_I \frac{\mathbf{C}\mathbf{b}_E}{\sqrt{\mathbf{b}_E^{\prime}\mathbf{P}\mathbf{b}_E}}}
#'
#' \strong{Implied economic weights:}
#' \deqn{\mathbf{w}_E = \frac{\sqrt{\lambda_E^2}}{\mathbf{b}_E^{\prime}\mathbf{P}\mathbf{b}_E} \mathbf{C}^{-1}\mathbf{P}\mathbf{b}_E}
#'
#' Uses \code{cpp_symmetric_solve} and \code{cpp_quadratic_form_sym} from
#' \code{math_primitives.cpp} for efficient matrix operations, and R's
#' \code{eigen()} for the eigendecomposition.
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 7.1.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#'
#' result <- esim(pmat, gmat)
#' print(result)
#' summary(result)
#' }
esim <- function(pmat, gmat, selection_intensity = 2.063, n_indices = 1L) {

  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  n_traits <- nrow(pmat)

  if (!isSymmetric(unname(pmat), tol = 1e-8))
    stop("pmat must be a symmetric matrix.")
  if (!isSymmetric(unname(gmat), tol = 1e-8))
    stop("gmat must be a symmetric matrix.")
  if (nrow(pmat) != nrow(gmat))
    stop("pmat and gmat must have the same dimensions.")
  if (n_traits < 2)
    stop("At least 2 traits are required for ESIM.")

  trait_names <- colnames(pmat)
  if (is.null(trait_names))
    trait_names <- paste0("Trait_", seq_len(n_traits))

  n_indices <- max(1L, as.integer(n_indices))

  # --------------------------------------------------------------------------
  # Step 1: Multi-trait heritability matrix H = P^{-1}C  (C = gmat)
  #   Uses cpp_symmetric_solve column-by-column (from math_primitives.cpp)
  #   H b = lambda^2 b  <=>  (P^{-1}G - lambda^2 I) b = 0
  # --------------------------------------------------------------------------
  P_inv_G <- .esim_solve_sym_multi(pmat, gmat)   # P^{-1}G via LDLT

  # --------------------------------------------------------------------------
  # Step 2: Eigendecomposition of P^{-1}G
  # --------------------------------------------------------------------------
  ev_result <- .leading_eigenvector(P_inv_G)

  lambda2  <- ev_result$value
  b_E      <- ev_result$vector
  all_vals <- ev_result$all_values

  # Collect up to n_indices leading vectors for summary
  all_pos   <- which(all_vals > 1e-8)
  n_avail   <- min(n_indices, length(all_pos))
  ev_full   <- eigen(P_inv_G, symmetric = FALSE)
  all_vecs  <- Re(ev_full$vectors)
  all_vals2 <- Re(ev_full$values)

  # --------------------------------------------------------------------------
  # Step 3: Compute metrics for primary index (b_E)  using C++ primitives
  # --------------------------------------------------------------------------
  metrics <- .eigen_index_metrics(b_E, pmat, gmat,
                                   lambda2 = lambda2,
                                   k_I    = selection_intensity)

  # --------------------------------------------------------------------------
  # Step 4: Implied economic weights
  #   w_E = (sqrt(lambda^2) / sigma_I^2) * C^{-1} P b_E
  #   = (|lambda| / bPb) * ginv(G) %*% P %*% b_E
  # --------------------------------------------------------------------------
  implied_w <- tryCatch({
    G_inv_Pb <- ginv(gmat) %*% (pmat %*% b_E)
    scale    <- sqrt(lambda2) / metrics$bPb
    as.numeric(scale * G_inv_Pb)
  }, error = function(e) {
    warning("Could not compute implied economic weights: ", e$message)
    rep(NA_real_, n_traits)
  })

  names(implied_w) <- trait_names

  # --------------------------------------------------------------------------
  # Step 5: Build summary data frame (one row per index requested)
  # --------------------------------------------------------------------------
  summary_rows <- vector("list", n_avail)
  pos_sorted   <- all_pos[order(all_vals2[all_pos], decreasing = TRUE)]

  for (j in seq_len(n_avail)) {
    idx_j  <- pos_sorted[j]
    b_j    <- Re(ev_full$vectors[, idx_j])
    b_j    <- b_j * sign(b_j[which.max(abs(b_j))])
    lam2_j <- all_vals2[idx_j]
    met_j  <- .eigen_index_metrics(b_j, pmat, gmat,
                                    lambda2 = lam2_j,
                                    k_I    = selection_intensity)
    b_row  <- setNames(round(b_j, 6),
                       paste0("b.", seq_len(n_traits)))
    summary_rows[[j]] <- data.frame(
      Index    = j,
      lambda2  = round(lam2_j, 6),
      hI2      = round(met_j$hI2, 6),
      rHI      = round(met_j$rHI, 6),
      sigma_I  = round(met_j$sigma_I, 6),
      Delta_G  = round(met_j$Delta_G, 6),
      as.data.frame(t(b_row)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  # --------------------------------------------------------------------------
  # Step 6: Assemble result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    b                   = setNames(b_E, trait_names),
    Delta_G             = setNames(metrics$Delta_G_vec, trait_names),
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    lambda2             = lambda2,
    implied_w           = implied_w,
    all_eigenvalues     = all_vals2,
    selection_intensity = selection_intensity,
    trait_names         = trait_names,
    pmat                = pmat,
    gmat                = gmat
  )

  class(result) <- c("esim", "eigen_index", "list")
  result
}

# ==============================================================================
# 7.2  RESIM - Restricted Eigen Selection Index
# ==============================================================================

#' Restricted Linear Phenotypic Eigen Selection Index (RESIM)
#'
#' @description
#' Extends ESIM by imposing null restrictions: genetic gains for \eqn{r} selected
#' traits are forced to zero while the index heritability for the remaining traits
#' is maximised.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param restricted_traits Integer vector of trait indices to restrict to zero
#'   genetic gain.  Example: \code{c(1, 3)} restricts traits 1 and 3.
#'   Alternatively supply a custom restriction matrix via \code{U_mat}.
#' @param U_mat Optional. Restriction matrix (n_traits x r) where each column
#'   defines one null restriction (\eqn{\mathbf{U}^{\prime}\mathbf{C}\mathbf{b} = 0}).
#'   Ignored if \code{restricted_traits} is provided.
#' @param selection_intensity Selection intensity constant (default: 2.063).
#'
#' @return Object of class \code{"resim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with b coefficients and key metrics.}
#'   \item{\code{b}}{Named numeric vector of RESIM coefficients.}
#'   \item{\code{Delta_G}}{Named vector of expected genetic gains per trait.}
#'   \item{\code{sigma_I}}{Index standard deviation.}
#'   \item{\code{hI2}}{Index heritability (leading eigenvalue of KP^(-1)C).}
#'   \item{\code{rHI}}{Index accuracy.}
#'   \item{\code{lambda2}}{Leading eigenvalue of the restricted eigenproblem.}
#'   \item{\code{K}}{Projection matrix used.}
#'   \item{\code{U_mat}}{Restriction matrix used.}
#'   \item{\code{restricted_traits}}{Integer vector of restricted trait indices.}
#'   \item{\code{implied_w}}{Implied economic weights.}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Projection matrix (Section 7.2):}
#' \deqn{\mathbf{K} = \mathbf{I}_t -
#'   \mathbf{P}^{-1}\mathbf{C}\mathbf{U}
#'   (\mathbf{U}^{\prime}\mathbf{C}\mathbf{P}^{-1}\mathbf{C}\mathbf{U})^{-1}
#'   \mathbf{U}^{\prime}\mathbf{C}}
#'
#' \strong{Restricted eigenproblem:}
#' \deqn{(\mathbf{K}\mathbf{P}^{-1}\mathbf{C} - \lambda_R^2 \mathbf{I}_t)\mathbf{b}_R = 0}
#'
#' \strong{Selection response and genetic gain:}
#' \deqn{R_R = k_I \sqrt{\mathbf{b}_R^{\prime}\mathbf{P}\mathbf{b}_R}}
#' \deqn{\mathbf{E}_R = k_I \frac{\mathbf{C}\mathbf{b}_R}{\sqrt{\mathbf{b}_R^{\prime}\mathbf{P}\mathbf{b}_R}}}
#'
#' \strong{Implied economic weights:}
#' \deqn{\mathbf{w}_R = \mathbf{C}^{-1}[\mathbf{P} + \mathbf{Q}_R^{\prime}\mathbf{C}]\mathbf{b}_R}
#' where \eqn{\mathbf{Q}_R = \mathbf{I} - \mathbf{K}}.
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 7.2.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#'
#' # Restrict traits 1 and 3 to zero genetic gain
#' result <- resim(pmat, gmat, restricted_traits = c(1, 3))
#' print(result)
#' summary(result)
#' }
resim <- function(pmat, gmat,
                  restricted_traits = NULL,
                  U_mat             = NULL,
                  selection_intensity = 2.063) {

  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  n_traits <- nrow(pmat)

  if (!isSymmetric(unname(pmat), tol = 1e-8))
    stop("pmat must be a symmetric matrix.")
  if (!isSymmetric(unname(gmat), tol = 1e-8))
    stop("gmat must be a symmetric matrix.")
  if (nrow(pmat) != nrow(gmat))
    stop("pmat and gmat must have the same dimensions.")

  trait_names <- colnames(pmat)
  if (is.null(trait_names))
    trait_names <- paste0("Trait_", seq_len(n_traits))

  # Build U matrix from restricted_traits or validate supplied U_mat
  if (!is.null(restricted_traits)) {
    if (!is.numeric(restricted_traits) ||
        length(restricted_traits) == 0L ||
        any(restricted_traits < 1) ||
        any(restricted_traits > n_traits)) {
      stop("restricted_traits must be valid trait indices between 1 and ", n_traits)
    }
    U_mat <- diag(n_traits)[, restricted_traits, drop = FALSE]
  } else if (is.null(U_mat)) {
    stop("Either 'restricted_traits' or 'U_mat' must be provided for RESIM.")
  } else {
    U_mat <- as.matrix(U_mat)
    if (nrow(U_mat) != n_traits)
      stop("U_mat must have nrow equal to n_traits (", n_traits, ")")
  }
  r_restrict <- ncol(U_mat)

  if (r_restrict >= n_traits)
    stop("Number of restrictions (", r_restrict,
         ") must be less than number of traits (", n_traits, ").")

  # --------------------------------------------------------------------------
  # Step 1: P^{-1}G and Projection matrix K
  #   K = I - P^{-1}G U (U'G P^{-1}G U)^{-1} U'G
  #   Using C++ cpp_symmetric_solve for each column of G and GU
  # --------------------------------------------------------------------------
  P_inv_G  <- .esim_solve_sym_multi(pmat, gmat)   # n_traits x n_traits
  GU       <- gmat %*% U_mat                       # G * U  (n_traits x r)
  P_inv_GU <- .esim_solve_sym_multi(pmat, GU)      # P^{-1}GU (n_traits x r)

  # Middle term: (U'G P^{-1}G U) = U'G (P^{-1}G U) = (GU)' (P^{-1}GU)
  mid  <- t(GU) %*% P_inv_GU                       # r x r
  # Solve mid^{-1} via ginv (handles near-singular cases)
  mid_inv <- ginv(mid)                              # r x r

  # Q_R = P^{-1}GU (mid_inv) U'G  (n_traits x n_traits)
  Q_R <- P_inv_GU %*% mid_inv %*% t(GU)
  K   <- diag(n_traits) - Q_R                      # projection matrix

  # --------------------------------------------------------------------------
  # Step 2: Restricted eigenproblem K P^{-1}G
  # --------------------------------------------------------------------------
  KPG <- K %*% P_inv_G

  ev_result <- .leading_eigenvector(KPG)
  lambda2   <- ev_result$value
  b_R       <- ev_result$vector

  # --------------------------------------------------------------------------
  # Step 3: Compute metrics via C++ primitives
  # --------------------------------------------------------------------------
  metrics <- .eigen_index_metrics(b_R, pmat, gmat,
                                   lambda2 = lambda2,
                                   k_I    = selection_intensity)

  # --------------------------------------------------------------------------
  # Step 4: Implied economic weights (Section 7.2)
  #   w_R = C^{-1}[P + Q_R' C] b_R
  # --------------------------------------------------------------------------
  implied_w <- tryCatch({
    inner   <- pmat + t(Q_R) %*% gmat
    G_inv   <- ginv(gmat)
    as.numeric(G_inv %*% (inner %*% b_R))
  }, error = function(e) {
    warning("Could not compute implied weights: ", e$message)
    rep(NA_real_, n_traits)
  })
  names(implied_w) <- trait_names

  # --------------------------------------------------------------------------
  # Step 5: Build summary data frame
  # --------------------------------------------------------------------------
  b_row      <- round(b_R, 6)
  b_cols     <- setNames(as.data.frame(t(b_row)),
                         paste0("b.", seq_len(n_traits)))
  summary_df <- data.frame(
    lambda2  = round(lambda2, 6),
    hI2      = round(metrics$hI2, 6),
    rHI      = round(metrics$rHI, 6),
    sigma_I  = round(metrics$sigma_I, 6),
    Delta_G  = round(metrics$Delta_G, 6),
    b_cols,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(summary_df) <- NULL

  # --------------------------------------------------------------------------
  # Step 6: Assemble result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    b                   = setNames(b_R, trait_names),
    Delta_G             = setNames(metrics$Delta_G_vec, trait_names),
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    lambda2             = lambda2,
    K                   = K,
    Q_R                 = Q_R,
    U_mat               = U_mat,
    restricted_traits   = restricted_traits,
    implied_w           = implied_w,
    selection_intensity = selection_intensity,
    trait_names         = trait_names,
    pmat                = pmat,
    gmat                = gmat
  )

  class(result) <- c("resim", "eigen_index", "list")
  result
}

# ==============================================================================
# 7.3  PPG-ESIM - Predetermined Proportional Gain Eigen Selection Index
# ==============================================================================

#' Predetermined Proportional Gain Eigen Selection Index (PPG-ESIM)
#'
#' @description
#' Extends ESIM by enforcing that genetic gains are proportional to a
#' user-specified vector \eqn{\mathbf{d}}: \eqn{\Delta\mathbf{G} \propto \mathbf{d}}.
#' A similarity transformation \eqn{\boldsymbol{\beta}_P = \mathbf{F}\mathbf{b}_P}
#' aligns the eigenvector with the breeder's desired direction.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param d Numeric vector of desired proportional gains (length n_traits).
#'   The \emph{ratios} among elements define target gain proportions.
#'   Direction (positive/negative) must reflect desired improvement direction
#'   (positive = increase, negative = decrease).
#' @param selection_intensity Selection intensity constant (default: 2.063).
#'
#' @return Object of class \code{"ppg_esim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with beta (transformed b), b (raw), hI2,
#'     rHI, sigma_I, Delta_G, and lambda2.}
#'   \item{\code{beta}}{Named numeric vector of post-transformation PPG-ESIM
#'     coefficients \eqn{\boldsymbol{\beta}_P = \mathbf{F}\mathbf{b}_P}.}
#'   \item{\code{b}}{Raw eigenvector b_P before similarity transform.}
#'   \item{\code{Delta_G}}{Named vector of expected genetic gains per trait.}
#'   \item{\code{sigma_I}}{Index standard deviation.}
#'   \item{\code{hI2}}{Index heritability.}
#'   \item{\code{rHI}}{Index accuracy.}
#'   \item{\code{lambda2}}{Leading eigenvalue of the PPG restricted eigenproblem.}
#'   \item{\code{F_mat}}{Diagonal similarity transform matrix F (diag(sign(d))).}
#'   \item{\code{K_P}}{PPG projection matrix (rank 1: projects onto d subspace).}
#'   \item{\code{D_M}}{Mallard matrix (t x t-1): orthogonal complement of d,
#'     used to construct the (t-1) restrictions.}
#'   \item{\code{desired_gains}}{Input proportional gains vector d.}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Restriction structure via the Mallard Matrix (Section 7.3):}
#'
#' The PPG-ESIM restricts the \eqn{(t-1)} directions **orthogonal** to \eqn{\mathbf{d}},
#' forcing the genetic gain vector to be collinear with \eqn{\mathbf{d}}.
#'
#' The Mallard matrix \eqn{\mathbf{D}_M} is \eqn{t \times (t-1)}: its columns span
#' the orthogonal complement of \eqn{\mathbf{d}}, obtained via QR decomposition of
#' \eqn{\mathbf{d}/\|\mathbf{d}\|}:
#' \deqn{\mathbf{Q}_{QR} = [\hat{d} \mid \mathbf{D}_M], \quad
#'       \text{QR}(\hat{d}) \to \mathbf{Q}_{QR} \in \mathbb{R}^{t \times t}}
#'
#' With \eqn{\boldsymbol{\Psi} = \mathbf{C}} (full-trait case, \eqn{\mathbf{U} = \mathbf{I}_t}):
#'
#' \strong{PPG projection matrix (\eqn{t-1} restrictions):}
#' \deqn{\mathbf{Q}_P =
#'   \mathbf{P}^{-1}\boldsymbol{\Psi}\mathbf{D}_M
#'   (\mathbf{D}_M^{\prime}\boldsymbol{\Psi}^{\prime}\mathbf{P}^{-1}\boldsymbol{\Psi}\mathbf{D}_M)^{-1}
#'   \mathbf{D}_M^{\prime}\boldsymbol{\Psi}^{\prime}}
#' \deqn{\mathbf{K}_P = \mathbf{I}_t - \mathbf{Q}_P \quad (\text{rank 1})}
#'
#' Because \eqn{\mathbf{K}_P} has rank 1 (projects onto the \eqn{\mathbf{d}} subspace),
#' \eqn{\mathbf{K}_P\mathbf{P}^{-1}\mathbf{C}} has exactly one positive eigenvalue and
#' its eigenvector produces \eqn{\Delta\mathbf{G} \propto \mathbf{d}}.
#'
#' \strong{PPG eigenproblem (rank-1 system):}
#' \deqn{(\mathbf{K}_P\mathbf{P}^{-1}\mathbf{C} - \lambda_P^2\mathbf{I}_t)\mathbf{b}_P = 0}
#'
#' \strong{Similarity transform:}
#' \deqn{\boldsymbol{\beta}_P = \mathbf{F}\mathbf{b}_P}
#' where \eqn{\mathbf{F} = \text{diag}(\text{sign}(\mathbf{d}))} aligns the
#' eigenvector sign with the breeder's intended improvement direction.
#'
#' \strong{Key response metrics:}
#' \deqn{R_P = k_I\sqrt{\boldsymbol{\beta}_P^{\prime}\mathbf{P}\boldsymbol{\beta}_P}}
#' \deqn{\mathbf{E}_P = k_I\frac{\mathbf{C}\boldsymbol{\beta}_P}{\sqrt{\boldsymbol{\beta}_P^{\prime}\mathbf{P}\boldsymbol{\beta}_P}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 7.3.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#' pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#'
#' # Desired proportional gains: increase all traits proportionally
#' d <- c(2, 1, 1, 1, 1, 1, 1)
#' result <- ppg_esim(pmat, gmat, d)
#' print(result)
#' summary(result)
#' }
ppg_esim <- function(pmat, gmat, d, selection_intensity = 2.063) {

  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  d    <- as.numeric(d)
  n_traits <- nrow(pmat)

  if (!isSymmetric(unname(pmat), tol = 1e-8))
    stop("pmat must be a symmetric matrix.")
  if (!isSymmetric(unname(gmat), tol = 1e-8))
    stop("gmat must be a symmetric matrix.")
  if (nrow(pmat) != nrow(gmat))
    stop("pmat and gmat must have the same dimensions.")
  if (length(d) != n_traits)
    stop("d must have the same length as the number of traits (", n_traits, ").")
  if (all(d == 0))
    stop("d must have at least one non-zero element.")

  trait_names <- colnames(pmat)
  if (is.null(trait_names))
    trait_names <- paste0("Trait_", seq_len(n_traits))

  # --------------------------------------------------------------------------
  # Step 1: Mallard matrix D_M - orthogonal complement of d
  #
  #   Goal: restrict the (t-1) directions ORTHOGONAL to d so that the
  #   remaining free direction forces Delta_G proportional to d.
  #
  #   Build via QR decomposition of the unit vector d/||d||:
  #     qr(d_unit) -> Q_full = [d_unit | D_M]
  #   D_M is t x (t-1) and its columns span orth(d).
  # --------------------------------------------------------------------------
  d_unit   <- d / sqrt(sum(d^2))                       # normalise d
  qr_d     <- qr(matrix(d_unit, ncol = 1))
  Q_full   <- qr.Q(qr_d, complete = TRUE)              # t x t orthogonal matrix
  D_M      <- Q_full[, -1L, drop = FALSE]              # t x (t-1): Mallard matrix

  # Psi = C * U  (U = I_t for full-trait PPG-ESIM, so Psi = G)
  Psi      <- gmat                                      # t x t

  # --------------------------------------------------------------------------
  # Step 2: Projection matrix K_P = I - Q_P  (rank t-1 restriction)
  #
  #   Q_P = P^{-1}Psi D_M (D_M'Psi'P^{-1}Psi D_M)^{-1} D_M'Psi'
  #
  #   Restricts the (t-1) dimensions orthogonal to d.
  #   K_P has rank 1: it projects onto the d subspace only.
  # --------------------------------------------------------------------------
  Psi_DM      <- Psi %*% D_M                              # t x (t-1)
  P_inv_PsiDM <- .esim_solve_sym_multi(pmat, Psi_DM)      # t x (t-1), via cpp_symmetric_solve

  # Middle matrix: D_M'Psi'P^{-1}Psi D_M  [(t-1) x (t-1)]
  mid    <- t(Psi_DM) %*% P_inv_PsiDM
  mid_inv <- tryCatch(
    solve(mid),
    error = function(e) {
      warning("Middle matrix in Q_P is singular; falling back to ginv(). ",
              "Check that d is not in the null space of P^{-1}G.", call. = FALSE)
      ginv(mid)
    }
  )

  # Q_P  [t x t, rank (t-1)]
  Q_P <- P_inv_PsiDM %*% mid_inv %*% t(Psi_DM)
  K_P <- diag(n_traits) - Q_P                            # rank-1 projection onto d

  # --------------------------------------------------------------------------
  # Step 3: PPG eigenproblem  (K_P P^{-1}G) b_P = lambda_P^2 b_P
  #
  #   Because K_P has rank 1 (projects onto the d subspace), we bypass the
  #   general eigendecomposition and solve analytically:
  #
  #     b_P  proportional to  K_P * P^{-1}G * d_unit   (any K_P*(non-null) gives the image)
  #     lambda_P^2 =  b_P' * P^{-1}G * b_P     (Rayleigh quotient; K_P*b_P = b_P)
  #
  #   This is numerically robust for rank-1 systems and avoids tolerance
  #   sensitivity in general eigendecomposition.
  # --------------------------------------------------------------------------
  P_inv_G    <- .esim_solve_sym_multi(pmat, gmat)     # t x t

  # Project P^{-1}G d_unit through K_P to find the single image direction
  v_raw  <- K_P %*% (P_inv_G %*% d_unit)             # t x 1
  v_norm <- sqrt(sum(v_raw^2))

  if (v_norm < 1e-12)
    stop("PPG-ESIM: K_P * P^{-1}G * d is numerically zero. ",
         "Check that d is not in the null space of the combined system.")

  b_P     <- as.numeric(v_raw) / v_norm               # normalised eigenvector

  # Eigenvalue: Rayleigh quotient of P^{-1}G at b_P
  # (valid because b_P lives in col(K_P) and K_P*b_P = b_P)
  lambda2 <- as.numeric(t(b_P) %*% P_inv_G %*% b_P)

  # --------------------------------------------------------------------------
  # Step 4: Similarity transformation beta_P = F b_P
  #
  #   The sign of b_P is arbitrary (eigenvector sign convention).
  #   We want Delta_G proportional to +d  (gains aligned with the breeder's desired direction).
  #   Check by computing the tentative gain direction G*b_P and comparing to d:
  #   if their dot product is negative, flip beta_P = -b_P.
  #
  #   This is the scalar version of the Mallard similarity transform F.
  # --------------------------------------------------------------------------
  tentative_gain <- as.numeric(gmat %*% b_P)   # G*b_P  (proportional to Delta_G)
  sign_corr      <- sign(sum(tentative_gain * d))  # +1 if Delta_G aligned with d
  if (sign_corr == 0L) sign_corr <- 1L
  beta_P <- b_P * sign_corr

  # Keep F_mat in output for reference (scalar +/-I in the rank-1 case)
  F_mat  <- sign_corr * diag(n_traits)

  # --------------------------------------------------------------------------
  # Step 5: Metrics computed on beta_P (the final coefficients after transform)
  # --------------------------------------------------------------------------
  metrics <- .eigen_index_metrics(beta_P, pmat, gmat,
                                   lambda2 = lambda2,
                                   k_I    = selection_intensity)

  # --------------------------------------------------------------------------
  # Step 6: Build summary data frame
  # --------------------------------------------------------------------------
  b_row    <- round(beta_P, 6)
  b_raw    <- round(b_P, 6)
  b_cols_beta <- setNames(as.data.frame(t(b_row)),  paste0("beta.", seq_len(n_traits)))
  b_cols_raw  <- setNames(as.data.frame(t(b_raw)),  paste0("b.",    seq_len(n_traits)))

  summary_df <- data.frame(
    lambda2  = round(lambda2, 6),
    hI2      = round(metrics$hI2, 6),
    rHI      = round(metrics$rHI, 6),
    sigma_I  = round(metrics$sigma_I, 6),
    Delta_G  = round(metrics$Delta_G, 6),
    b_cols_beta,
    b_cols_raw,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(summary_df) <- NULL

  # --------------------------------------------------------------------------
  # Step 7: Assemble result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    beta                = setNames(beta_P, trait_names),
    b                   = setNames(b_P, trait_names),
    Delta_G             = setNames(metrics$Delta_G_vec, trait_names),
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    lambda2             = lambda2,
    F_mat               = F_mat,
    K_P                 = K_P,
    Q_P                 = Q_P,
    D_M                 = D_M,        # Mallard matrix (t x t-1): orth complement of d
    desired_gains       = setNames(d, trait_names),
    selection_intensity = selection_intensity,
    trait_names         = trait_names,
    pmat                = pmat,
    gmat                = gmat
  )

  class(result) <- c("ppg_esim", "eigen_index", "list")
  result
}

# ==============================================================================
# S3 PRINT METHODS
# ==============================================================================

#' Print method for ESIM
#' @param x Object of class 'esim'
#' @param ... Additional arguments (unused)
#' @export
print.esim <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("LINEAR PHENOTYPIC EIGEN SELECTION INDEX (ESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 7, Section 7.1\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("LEADING EIGENVALUE (Maximised Index Heritability)\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_E^2  (h^2_I_E):  %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy r_HI:          %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev sigma_I:  %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response R_E: %.6f\n", x$sigma_I * x$selection_intensity))

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("ESIM INDEX COEFFICIENTS (b_E)\n")
  cat("-------------------------------------------------------------\n")
  b_df <- data.frame(
    Trait       = x$trait_names,
    Coefficient = round(as.numeric(x$b), 6),
    stringsAsFactors = FALSE
  )
  print(b_df, row.names = FALSE)

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_E)\n")
  cat("-------------------------------------------------------------\n")
  eg_df <- data.frame(
    Trait   = x$trait_names,
    Delta_G = round(as.numeric(x$Delta_G), 6),
    stringsAsFactors = FALSE
  )
  print(eg_df, row.names = FALSE)

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("IMPLIED ECONOMIC WEIGHTS (w_E)\n")
  cat("-------------------------------------------------------------\n")
  iw_df <- data.frame(
    Trait    = x$trait_names,
    Implied_w = round(as.numeric(x$implied_w), 6),
    stringsAsFactors = FALSE
  )
  print(iw_df, row.names = FALSE)

  if (length(x$all_eigenvalues) > 1) {
    cat("\n")
    cat("-------------------------------------------------------------\n")
    cat("ALL REAL EIGENVALUES OF P^{-1}G\n")
    cat("-------------------------------------------------------------\n")
    real_vals <- Re(x$all_eigenvalues)
    cat(round(sort(real_vals, decreasing = TRUE), 6), "\n")
  }

  cat("\n")
  invisible(x)
}

#' Summary method for ESIM
#' @param object Object of class 'esim'
#' @param ... Additional arguments (unused)
#' @export
summary.esim <- function(object, ...) {
  print(object, ...)

  cat("\n")
  cat("==============================================================\n")
  cat("INDEX SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")

  cat("Notes:\n")
  cat("  - lambda^2 = eigenvalue = maximised index heritability h^2_I_E\n")
  cat("  - b_E are optimal without requiring pre-specified economic weights\n")
  cat("  - Implied_w reconstructs equivalent economic weights post-hoc\n\n")

  invisible(object)
}


#' Print method for RESIM
#' @param x Object of class 'resim'
#' @param ... Additional arguments (unused)
#' @export
print.resim <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("RESTRICTED EIGEN SELECTION INDEX (RESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 7, Section 7.2\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n")

  if (!is.null(x$restricted_traits)) {
    cat("Restricted traits (zero gain):",
        paste(x$trait_names[x$restricted_traits], collapse = ", "), "\n")
  }
  cat("Number of restrictions:   ", ncol(x$U_mat), "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("RESIM LEADING EIGENVALUE\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_R^2  (h^2_I_R):  %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy r_HI:          %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev sigma_I:  %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response R_R: %.6f\n", x$sigma_I * x$selection_intensity))

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("RESIM INDEX COEFFICIENTS (b_R)\n")
  cat("-------------------------------------------------------------\n")
  b_df <- data.frame(
    Trait       = x$trait_names,
    Coefficient = round(as.numeric(x$b), 6),
    stringsAsFactors = FALSE
  )
  print(b_df, row.names = FALSE)

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_R)\n")
  cat("-------------------------------------------------------------\n")
  eg_df <- data.frame(
    Trait      = x$trait_names,
    Delta_G    = round(as.numeric(x$Delta_G), 6),
    Restricted = if (!is.null(x$restricted_traits)) {
      seq_along(x$trait_names) %in% x$restricted_traits
    } else {
      rep(FALSE, length(x$trait_names))
    },
    stringsAsFactors = FALSE
  )
  print(eg_df, row.names = FALSE)

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("IMPLIED ECONOMIC WEIGHTS (w_R)\n")
  cat("-------------------------------------------------------------\n")
  iw_df <- data.frame(
    Trait     = x$trait_names,
    Implied_w = round(as.numeric(x$implied_w), 6),
    stringsAsFactors = FALSE
  )
  print(iw_df, row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' Summary method for RESIM
#' @param object Object of class 'resim'
#' @param ... Additional arguments (unused)
#' @export
summary.resim <- function(object, ...) {
  print(object, ...)

  cat("\n")
  cat("==============================================================\n")
  cat("INDEX SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")

  cat("Notes:\n")
  cat("  - K is the projection matrix I - P^{-1}GU(U'GP^{-1}GU)^{-1}U'G\n")
  cat("  - lambda_R^2 is the leading eigenvalue of K*P^{-1}*G\n")
  cat("  - Restricted traits have |Delta_G| <= numerical tolerance\n\n")

  invisible(object)
}


#' Print method for PPG-ESIM
#' @param x Object of class 'ppg_esim'
#' @param ... Additional arguments (unused)
#' @export
print.ppg_esim <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("PREDETERMINED PROPORTIONAL GAIN EIGEN SELECTION INDEX (PPG-ESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 7, Section 7.3\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("PPG-ESIM LEADING EIGENVALUE\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_P^2  (h^2_I_P):  %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy r_HI:          %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev sigma_I:  %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response R_P: %.6f\n", x$sigma_I * x$selection_intensity))

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("PPG-ESIM COEFFICIENTS AFTER SIMILARITY TRANSFORM (beta_P = F b_P)\n")
  cat("-------------------------------------------------------------\n")
  b_df <- data.frame(
    Trait       = x$trait_names,
    beta_P      = round(as.numeric(x$beta), 6),
    b_P_raw     = round(as.numeric(x$b), 6),
    stringsAsFactors = FALSE
  )
  print(b_df, row.names = FALSE)

  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("DESIRED vs. ACHIEVED PROPORTIONAL GAINS\n")
  cat("-------------------------------------------------------------\n")
  Delta_G_vals <- as.numeric(x$Delta_G)
  desired_vals <- as.numeric(x$desired_gains)

  # Compute proportionality ratios
  ratios <- ifelse(abs(desired_vals) > 1e-10,
                   Delta_G_vals / desired_vals, NA_real_)
  avg_ratio <- mean(ratios, na.rm = TRUE)

  pg_df <- data.frame(
    Trait       = x$trait_names,
    Desired_d   = round(desired_vals, 6),
    Achieved_E  = round(Delta_G_vals, 6),
    Ratio       = round(ratios, 6),
    stringsAsFactors = FALSE
  )
  print(pg_df, row.names = FALSE)

  cv_ratio <- if (!is.na(avg_ratio) && avg_ratio != 0) {
    sd(ratios, na.rm = TRUE) / abs(avg_ratio)
  } else NA_real_

  cat(sprintf("\nMean proportionality scale (phi): %.6f\n", avg_ratio))
  cat("The Mallard-matrix construction constrains gains to be parallel to d.\n")
  if (!is.na(cv_ratio) && cv_ratio < 0.02) {
    cat(sprintf("[OK]  Proportional gains achieved (CV of ratios: %.4f%%)\n",
                cv_ratio * 100))
  } else if (!is.na(cv_ratio)) {
    cat(sprintf("[!]   Residual proportionality deviation (CV: %.4f%%).\n",
                cv_ratio * 100))
    cat("      Check matrix conditioning or whether d conflicts with\n")
    cat("      the genetic covariance structure.\n")
  }

  cat("\n")
  invisible(x)
}

#' Summary method for PPG-ESIM
#' @param object Object of class 'ppg_esim'
#' @param ... Additional arguments (unused)
#' @export
summary.ppg_esim <- function(object, ...) {
  print(object, ...)

  cat("\n")
  cat("==============================================================\n")
  cat("INDEX SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")

  cat("Notes:\n")
  cat("  - D_M (Mallard matrix, t x t-1): orthogonal complement of d\n")
  cat("  - Q_P = P^{-1}*Psi*D_M * (D_M'*Psi'*P^{-1}*Psi*D_M)^{-1} * D_M'*Psi'\n")
  cat("  - K_P = I - Q_P  (rank 1: projects onto d direction only)\n")
  cat("  - beta_P = F b_P  (similarity transform for sign alignment with d)\n")
  cat("  - Proportional gains: ratios Achieved/Desired should be ~constant\n\n")

  invisible(object)
}