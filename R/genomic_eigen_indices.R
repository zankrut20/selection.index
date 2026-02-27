#' Linear Molecular and Genomic Eigen Selection Index Methods (Chapter 8)
#' @name genomic_eigen_indices
#'
#' @description
#' Implements the Linear Molecular and Genomic Eigen Selection Index methods from Chapter 8.
#' These methods extend eigen-based selection to genomic/molecular data, maximizing
#' the accuracy squared (rho_HI^2) through eigenvalue problems without requiring
#' pre-specified economic weights.
#'
#' Methods included:
#' - MESIM    : Molecular Eigen Selection Index Method (Section 8.1)
#' - GESIM    : Linear Genomic Eigen Selection Index Method (Section 8.2)
#' - GW-ESIM  : Genome-Wide Linear Eigen Selection Index Method (Section 8.3)
#' - RGESIM   : Restricted Linear Genomic Eigen Selection Index Method (Section 8.4)
#' - PPG-GESIM: Predetermined Proportional Gain Genomic Eigen Selection Index (Section 8.5)
#'
#' All implementations use C++ primitives (math_primitives.cpp) for quadratic forms
#' and symmetric solves, while eigendecompositions use R's eigen() for correctness
#' and compatibility with the existing package architecture.
#'
#' @section Mathematical Foundation:
#'
#' Like the phenotypic ESIM (Chapter 7), these genomic eigen methods maximize
#' the squared accuracy between the index and net genetic merit, but incorporate
#' molecular markers, GEBVs, or genome-wide marker scores.
#'
#' The general approach solves a generalized eigenproblem to find optimal index
#' coefficients that maximize heritability without requiring economic weights.
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant
#' Breeding. Springer International Publishing. Chapter 8.
#'
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom MASS ginv
NULL

# ==============================================================================
# LOCAL HELPERS
# ==============================================================================

#' Solve symmetric linear system for multiple right-hand sides
#' @keywords internal
#' @noRd
.gesim_solve_sym_multi <- function(A, B) {
  B <- as.matrix(B)
  out <- matrix(0, nrow = nrow(A), ncol = ncol(B))
  for (j in seq_len(ncol(B))) {
    out[, j] <- cpp_symmetric_solve(A, B[, j])
  }
  out
}

#' Compute standard metrics for a genomic eigen-based index
#'
#' @description
#' Computes all standard metrics using C++ math primitives.
#' For genomic eigen indices, the eigenvalue lambda^2 IS the index heritability.
#'
#' @param b  Index coefficient vector (eigenvector or transformed eigenvector)
#' @param Phi  Combined phenotypic/marker variance-covariance matrix
#' @param A  Combined genotypic/marker variance-covariance matrix
#' @param lambda2  Eigenvalue associated with b (used for h^2_I when exact)
#' @param k_I  Selection intensity constant (default 2.063)
#' @keywords internal
#' @noRd
.genomic_eigen_index_metrics <- function(b, Phi, A, lambda2 = NULL, k_I = 2.063) {
  b <- as.numeric(b)

  # --- Quadratic forms via C++ primitives ---
  bPhibb <- cpp_quadratic_form_sym(b, Phi) # b'Phi*b
  bAb <- cpp_quadratic_form_sym(b, A) # b'A*b

  # Ensure positive quadratic form (eigenvector sign is arbitrary)
  if (bPhibb < 0) {
    b <- -b
    bPhibb <- -bPhibb
  }

  sigma_I <- if (bPhibb > 0) sqrt(bPhibb) else NA_real_

  # Selection response: R = k_I * sigma_I
  R <- if (!is.na(sigma_I)) k_I * sigma_I else NA_real_

  # Expected genetic gain per trait: E = (k_I / sigma_I) * A*b
  E_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    as.vector(k_I * (A %*% b) / sigma_I)
  } else {
    rep(NA_real_, nrow(A))
  }

  # Index heritability: use eigenvalue if available (exact), else b'Ab / b'Phib
  hI2 <- if (!is.null(lambda2)) {
    as.numeric(lambda2)
  } else if (!is.na(bPhibb) && bPhibb > 0) {
    bAb / bPhibb
  } else {
    NA_real_
  }

  # Accuracy: r_HI = sqrt(h^2_I)
  rHI <- if (!is.na(hI2) && hI2 >= 0) sqrt(hI2) else NA_real_

  list(
    b          = b, # Return corrected eigenvector
    bPhibb     = bPhibb,
    bAb        = bAb,
    sigma_I    = sigma_I,
    R          = R,
    E_vec      = E_vec,
    hI2        = hI2,
    rHI        = rHI
  )
}

#' Select the leading eigenvector from a real-eigendecomposition
#'
#' @description
#' Returns the eigenvector paired with the *largest real* eigenvalue that is
#' positive (negative eigenvalues indicate numerical noise or rank deficiency).
#' Normalises the sign so that the largest-magnitude element is positive.
#'
#' @param mat Square matrix to decompose
#' @param tol Eigenvalue tolerance (default 1e-8)
#' @return List: \code{vector}, \code{value}, \code{all_values}
#' @keywords internal
#' @noRd
.gesim_leading_eigenvector <- function(mat, tol = 1e-8) {
  ev <- eigen(mat, symmetric = FALSE)
  vals <- Re(ev$values) # work with real parts
  vecs <- Re(ev$vectors)

  # Keep only positive eigenvalues
  pos <- which(vals > tol)
  if (length(pos) == 0) {
    stop(
      "No positive eigenvalues found. ",
      "Check that matrices are valid variance-covariance matrices."
    )
  }

  # Leading (largest positive) eigenvalue
  idx <- pos[which.max(vals[pos])]
  bvec <- vecs[, idx]

  # Canonical sign: make the largest-magnitude element positive
  bvec <- bvec * sign(bvec[which.max(abs(bvec))])

  list(vector = bvec, value = vals[idx], all_values = vals)
}

# ==============================================================================
# 8.1  MESIM - Molecular Eigen Selection Index Method
# ==============================================================================

#' Molecular Eigen Selection Index Method (MESIM)
#'
#' @description
#' Implements the MESIM by maximising the squared accuracy through the
#' generalised eigenproblem combining phenotypic data with molecular marker scores.
#' Unlike Smith-Hazel LPSI, **no economic weights are required**.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param S_M Covariance between phenotypes and marker scores (n_traits x n_traits).
#'   This represents Cov(y, s) where s are marker scores. Used in the phenotypic
#'   variance matrix T_M.
#' @param S_Mg Covariance between genetic values and marker scores (n_traits x n_traits).
#'   This represents Cov(g, s). Used in the genetic variance matrix Psi_M.
#'   If NULL (default), uses S_M, which is appropriate when assuming
#'   Cov(e, s) ~= 0 (errors uncorrelated with markers), so Cov(y,s) ~= Cov(g,s).
#' @param S_var Variance-covariance matrix of marker scores (n_traits x n_traits).
#'   Represents Var(s). If NULL (default), uses S_M for backward compatibility,
#'   but providing the actual Var(s) is more statistically rigorous.
#' @param selection_intensity Selection intensity constant \eqn{k_I}
#'   (default: 2.063 for 10\% selection).
#'
#' @return Object of class \code{"mesim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with coefficients and metrics.}
#'   \item{\code{b_y}}{Coefficients for phenotypic data.}
#'   \item{\code{b_s}}{Coefficients for marker scores.}
#'   \item{\code{b_combined}}{Combined coefficient vector.}
#'   \item{\code{E_M}}{Expected genetic gains per trait.}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{hI2}}{Index heritability (= leading eigenvalue).}
#'   \item{\code{rHI}}{Accuracy \eqn{r_{HI}}.}
#'   \item{\code{R_M}}{Selection response.}
#'   \item{\code{lambda2}}{Leading eigenvalue (maximised index heritability).}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Eigenproblem (Section 8.1):}
#' \deqn{(\mathbf{T}_M^{-1}\mathbf{\Psi}_M - \lambda_M^2 \mathbf{I}_{2t})\boldsymbol{\beta}_M = 0}
#'
#' where:
#' \deqn{\mathbf{T}_M = \begin{bmatrix} \mathbf{P} & \mathrm{Cov}(\mathbf{y},\mathbf{s}) \\ \mathrm{Cov}(\mathbf{s},\mathbf{y}) & \mathrm{Var}(\mathbf{s}) \end{bmatrix}}
#' \deqn{\mathbf{\Psi}_M = \begin{bmatrix} \mathbf{C} & \mathrm{Cov}(\mathbf{g},\mathbf{s}) \\ \mathrm{Cov}(\mathbf{s},\mathbf{g}) & \mathrm{Var}(\mathbf{s}) \end{bmatrix}}
#'
#' \strong{Theoretical distinction:}
#' \itemize{
#'   \item T_M uses phenotypic covariances: Cov(y, s) provided via \code{S_M}
#'   \item Psi_M uses genetic covariances: Cov(g, s) provided via \code{S_Mg}
#'   \item Since y = g + e, if Cov(e, s) ~= 0, then Cov(y, s) ~= Cov(g, s)
#'   \item Chapter 8.1 assumes marker scores are pure genetic predictors,
#'         so S_M can be used for both (default behavior when S_Mg = NULL)
#' }
#'
#' The solution \eqn{\lambda_M^2} (largest eigenvalue) equals the maximum
#' achievable index heritability.
#'
#' \strong{Key metrics:}
#' \deqn{R_M = k_I \sqrt{\boldsymbol{\beta}_{M}^{\prime}\mathbf{T}_M\boldsymbol{\beta}_{M}}}
#' \deqn{\mathbf{E}_M = k_I \frac{\mathbf{\Psi}_M\boldsymbol{\beta}_{M}}{\sqrt{\boldsymbol{\beta}_{M}^{\prime}\mathbf{T}_M\boldsymbol{\beta}_{M}}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 8.1.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate marker score matrices (in practice, compute from data)
#' S_M <- gmat * 0.7 # Cov(y, s) - phenotype-marker covariance
#' S_Mg <- gmat * 0.65 # Cov(g, s) - genetic-marker covariance
#' S_var <- gmat * 0.8 # Var(s) - marker score variance
#'
#' # Most rigorous: Provide all three covariance matrices
#' result <- mesim(pmat, gmat, S_M, S_Mg = S_Mg, S_var = S_var)
#' print(result)
#'
#' # Standard usage: Cov(g,s) defaults to Cov(y,s) when errors uncorrelated
#' result_standard <- mesim(pmat, gmat, S_M, S_var = S_var)
#'
#' # Backward compatible: Chapter 8.1 simplified notation
#' result_simple <- mesim(pmat, gmat, S_M)
#' }
mesim <- function(pmat, gmat, S_M, S_Mg = NULL, S_var = NULL, selection_intensity = 2.063) {
  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  S_M <- as.matrix(S_M)
  n_traits <- nrow(pmat)

  # If S_Mg not provided, use S_M (assumes Cov(e,s) ~= 0, so Cov(y,s) ~= Cov(g,s))
  if (is.null(S_Mg)) {
    S_Mg <- S_M
  } else {
    S_Mg <- as.matrix(S_Mg)
  }

  # If S_var not provided, use S_M (backward compatibility with Chapter 8.1)
  if (is.null(S_var)) {
    S_var <- S_M
  } else {
    S_var <- as.matrix(S_var)
  }

  if (!isSymmetric(unname(pmat), tol = 1e-8)) {
    stop("pmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(gmat), tol = 1e-8)) {
    stop("gmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(S_M), tol = 1e-8)) {
    stop("S_M must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(S_Mg), tol = 1e-8)) {
    stop("S_Mg must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(S_var), tol = 1e-8)) {
    stop("S_var must be a symmetric matrix.")
  }
  if (nrow(pmat) != nrow(gmat) || nrow(pmat) != nrow(S_M) ||
    nrow(pmat) != nrow(S_Mg) || nrow(pmat) != nrow(S_var)) {
    stop("All matrices must have the same dimensions.")
  }
  if (n_traits < 2) {
    stop("At least 2 traits are required for MESIM.")
  }

  trait_names <- colnames(pmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  # --------------------------------------------------------------------------
  # Step 1: Construct combined matrices T_M and Psi_M
  # --------------------------------------------------------------------------
  # THEORETICAL STRUCTURE (Section 8.1):
  #
  # T_M (Phenotypic variance matrix):
  #   [Var(y)      Cov(y,s)]   = [P      S_M  ]
  #   [Cov(s,y)    Var(s)  ]     [S_M'   S_var]
  #
  # Psi_M (Genetic variance matrix):
  #   [Var_G(y)    Cov(g,s)]   = [C      S_Mg ]
  #   [Cov(s,g)    Var(s)  ]     [S_Mg'  S_var]
  #
  # KEY DISTINCTION:
  # - T_M uses Cov(y,s) because it describes phenotypic covariance structure
  # - Psi_M uses Cov(g,s) because it describes genetic covariance structure
  # - When marker scores are pure genetic predictors and errors are uncorrelated
  #   with markers (Cov(e,s) ~= 0), then Cov(y,s) ~= Cov(g,s), so S_Mg ~= S_M
  #   (this is the Chapter 8.1 textbook assumption)

  T_M <- rbind(
    cbind(pmat, S_M),
    cbind(S_M, S_var)
  )

  Psi_M <- rbind(
    cbind(gmat, S_Mg), # Use Cov(g,s), not Cov(y,s)
    cbind(S_Mg, S_var)
  )

  # --------------------------------------------------------------------------
  # Step 2: Solve eigenproblem T_M^{-1} Psi_M
  # --------------------------------------------------------------------------
  T_M_inv_Psi_M <- .gesim_solve_sym_multi(T_M, Psi_M)

  ev_result <- .gesim_leading_eigenvector(T_M_inv_Psi_M)
  lambda2 <- ev_result$value
  b_M <- ev_result$vector

  # --------------------------------------------------------------------------
  # Step 3: Compute metrics using combined matrices
  # --------------------------------------------------------------------------
  metrics <- .genomic_eigen_index_metrics(b_M, T_M, Psi_M,
    lambda2 = lambda2,
    k_I = selection_intensity
  )

  # Extract corrected eigenvector (sign-corrected for positive quadratic form)
  b_M <- metrics$b

  # Split into phenotype and marker score coefficients
  b_y <- b_M[1:n_traits]
  b_s <- b_M[(n_traits + 1):(2 * n_traits)]

  # Expected gains are first n_traits elements of E_vec
  E_M <- metrics$E_vec[1:n_traits]
  names(E_M) <- trait_names

  # --------------------------------------------------------------------------
  # Step 4: Build summary data frame
  # --------------------------------------------------------------------------
  b_y_vec <- round(b_y, 6)
  b_s_vec <- round(b_s, 6)

  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_len(n_traits))

  b_s_df <- as.data.frame(matrix(b_s_vec, nrow = 1))
  colnames(b_s_df) <- paste0("b_s.", seq_len(n_traits))

  summary_df <- data.frame(
    b_y_df,
    b_s_df,
    hI2 = round(metrics$hI2, 6),
    rHI = round(metrics$rHI, 6),
    sigma_I = round(metrics$sigma_I, 6),
    R_M = round(metrics$R, 6),
    lambda2 = round(lambda2, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # --------------------------------------------------------------------------
  # Step 5: Return result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    b_y                 = setNames(b_y, trait_names),
    b_s                 = setNames(b_s, trait_names),
    b_combined          = b_M,
    E_M                 = E_M,
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    R_M                 = metrics$R,
    lambda2             = lambda2,
    selection_intensity = selection_intensity,
    T_M                 = T_M,
    Psi_M               = Psi_M,
    trait_names         = trait_names
  )

  class(result) <- c("mesim", "genomic_eigen_index", "list")
  result
}

# ==============================================================================
# 8.2  GESIM - Linear Genomic Eigen Selection Index Method
# ==============================================================================

#' Linear Genomic Eigen Selection Index Method (GESIM)
#'
#' @description
#' Implements the GESIM by maximising the squared accuracy through the
#' generalised eigenproblem combining phenotypic data with GEBVs (Genomic
#' Estimated Breeding Values). No economic weights are required.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param Gamma Covariance between phenotypes and GEBVs (n_traits x n_traits).
#'   This represents Cov(y, gamma) where gamma are GEBVs.
#' @param selection_intensity Selection intensity constant \eqn{k_I}
#'   (default: 2.063 for 10\% selection).
#'
#' @return Object of class \code{"gesim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with coefficients and metrics.}
#'   \item{\code{b_y}}{Coefficients for phenotypic data.}
#'   \item{\code{b_gamma}}{Coefficients for GEBVs.}
#'   \item{\code{b_combined}}{Combined coefficient vector.}
#'   \item{\code{E_G}}{Expected genetic gains per trait.}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{hI2}}{Index heritability (= leading eigenvalue).}
#'   \item{\code{rHI}}{Accuracy \eqn{r_{HI}}.}
#'   \item{\code{R_G}}{Selection response.}
#'   \item{\code{lambda2}}{Leading eigenvalue.}
#'   \item{\code{implied_w}}{Implied economic weights.}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Eigenproblem (Section 8.2):}
#' \deqn{(\mathbf{\Phi}^{-1}\mathbf{A} - \lambda_G^2 \mathbf{I}_{2t})\boldsymbol{\beta}_G = 0}
#'
#' where:
#' \deqn{\mathbf{\Phi} = \begin{bmatrix} \mathbf{P} & \mathbf{\Gamma} \\ \mathbf{\Gamma} & \mathbf{\Gamma} \end{bmatrix}}
#' \deqn{\mathbf{A} = \begin{bmatrix} \mathbf{C} & \mathbf{\Gamma} \\ \mathbf{\Gamma} & \mathbf{\Gamma} \end{bmatrix}}
#'
#' \strong{Implied economic weights:}
#' \deqn{\mathbf{w}_G = \mathbf{A}^{-1}\mathbf{\Phi}\boldsymbol{\beta}}
#'
#' \strong{Selection response:}
#' \deqn{R_G = k_I \sqrt{\boldsymbol{\beta}_G^{\prime}\mathbf{\Phi}\boldsymbol{\beta}_G}}
#'
#' \strong{Expected genetic gain per trait:}
#' \deqn{\mathbf{E}_G = k_I \frac{\mathbf{A}\boldsymbol{\beta}_G}{\sqrt{\boldsymbol{\beta}_G^{\prime}\mathbf{\Phi}\boldsymbol{\beta}_G}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 8.2.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate GEBV covariance (in practice, compute from genomic predictions)
#' Gamma <- gmat * 0.8 # Assume 80% GEBV-phenotype covariance
#'
#' result <- gesim(pmat, gmat, Gamma)
#' print(result)
#' }
gesim <- function(pmat, gmat, Gamma, selection_intensity = 2.063) {
  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  Gamma <- as.matrix(Gamma)
  n_traits <- nrow(pmat)

  if (!isSymmetric(unname(pmat), tol = 1e-8)) {
    stop("pmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(gmat), tol = 1e-8)) {
    stop("gmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(Gamma), tol = 1e-8)) {
    stop("Gamma must be a symmetric matrix.")
  }
  if (nrow(pmat) != nrow(gmat) || nrow(pmat) != nrow(Gamma)) {
    stop("All matrices must have the same dimensions.")
  }
  if (n_traits < 2) {
    stop("At least 2 traits are required for GESIM.")
  }

  trait_names <- colnames(pmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  # --------------------------------------------------------------------------
  # Step 1: Construct combined matrices Phi and A
  # --------------------------------------------------------------------------
  # Phi = [P      Gamma]
  #       [Gamma  Gamma]
  Phi <- rbind(
    cbind(pmat, Gamma),
    cbind(Gamma, Gamma)
  )

  # A = [C      Gamma]
  #     [Gamma  Gamma]
  A <- rbind(
    cbind(gmat, Gamma),
    cbind(Gamma, Gamma)
  )

  # --------------------------------------------------------------------------
  # Step 2: Solve eigenproblem Phi^{-1} A
  # --------------------------------------------------------------------------
  Phi_inv_A <- .gesim_solve_sym_multi(Phi, A)

  ev_result <- .gesim_leading_eigenvector(Phi_inv_A)
  lambda2 <- ev_result$value
  b_G <- ev_result$vector

  # --------------------------------------------------------------------------
  # Step 3: Compute metrics
  # --------------------------------------------------------------------------
  metrics <- .genomic_eigen_index_metrics(b_G, Phi, A,
    lambda2 = lambda2,
    k_I = selection_intensity
  )

  # Extract corrected eigenvector (sign-corrected for positive quadratic form)
  b_G <- metrics$b

  # Split into phenotype and GEBV coefficients
  b_y <- b_G[1:n_traits]
  b_gamma <- b_G[(n_traits + 1):(2 * n_traits)]

  # Expected gains are first n_traits elements
  E_G <- metrics$E_vec[1:n_traits]
  names(E_G) <- trait_names

  # --------------------------------------------------------------------------
  # Step 4: Implied economic weights: w_G = A^{-1} Phi beta
  # --------------------------------------------------------------------------
  implied_w <- tryCatch(
    {
      A_inv_Phi_b <- ginv(A) %*% (Phi %*% b_G)
      as.numeric(A_inv_Phi_b[1:n_traits])
    },
    error = function(e) {
      warning("Could not compute implied economic weights: ", e$message)
      rep(NA_real_, n_traits)
    }
  )
  names(implied_w) <- trait_names

  # --------------------------------------------------------------------------
  # Step 5: Build summary data frame
  # --------------------------------------------------------------------------
  b_y_vec <- round(b_y, 6)
  b_gamma_vec <- round(b_gamma, 6)

  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_len(n_traits))

  b_gamma_df <- as.data.frame(matrix(b_gamma_vec, nrow = 1))
  colnames(b_gamma_df) <- paste0("b_gamma.", seq_len(n_traits))

  summary_df <- data.frame(
    b_y_df,
    b_gamma_df,
    hI2 = round(metrics$hI2, 6),
    rHI = round(metrics$rHI, 6),
    sigma_I = round(metrics$sigma_I, 6),
    R_G = round(metrics$R, 6),
    lambda2 = round(lambda2, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # --------------------------------------------------------------------------
  # Step 6: Return result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    b_y                 = setNames(b_y, trait_names),
    b_gamma             = setNames(b_gamma, trait_names),
    b_combined          = b_G,
    E_G                 = E_G,
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    R_G                 = metrics$R,
    lambda2             = lambda2,
    implied_w           = implied_w,
    selection_intensity = selection_intensity,
    Phi                 = Phi,
    A                   = A,
    trait_names         = trait_names
  )

  class(result) <- c("gesim", "genomic_eigen_index", "list")
  result
}

# ==============================================================================
# 8.3  GW-ESIM - Genome-Wide Linear Eigen Selection Index Method
# ==============================================================================

#' Genome-Wide Linear Eigen Selection Index Method (GW-ESIM)
#'
#' @description
#' Implements the GW-ESIM by incorporating genome-wide marker effects directly
#' into the eigen selection index framework. Uses N marker scores alongside
#' phenotypic data.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param G_M Covariance between phenotypes and marker scores (n_traits x N_markers).
#' @param M Variance-covariance matrix of marker scores (N_markers x N_markers).
#' @param selection_intensity Selection intensity constant \eqn{k_I}
#'   (default: 2.063 for 10\% selection).
#'
#' @return Object of class \code{"gw_esim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with key metrics.}
#'   \item{\code{b_y}}{Coefficients for phenotypic data.}
#'   \item{\code{b_m}}{Coefficients for marker scores.}
#'   \item{\code{b_combined}}{Combined coefficient vector.}
#'   \item{\code{E_W}}{Expected genetic gains per trait.}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{hI2}}{Index heritability (= leading eigenvalue).}
#'   \item{\code{rHI}}{Accuracy.}
#'   \item{\code{R_W}}{Selection response.}
#'   \item{\code{lambda2}}{Leading eigenvalue.}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Eigenproblem (Section 8.3):}
#' \deqn{(\mathbf{Q}^{-1}\mathbf{X} - \lambda_W^2 \mathbf{I}_{(t+N)})\boldsymbol{\beta}_W = 0}
#'
#' where:
#' \deqn{\mathbf{Q} = \begin{bmatrix} \mathbf{P} & \mathbf{G}_M \\ \mathbf{G}_M^{\prime} & \mathbf{M} \end{bmatrix}}
#' \deqn{\mathbf{X} = \begin{bmatrix} \mathbf{C} & \mathbf{G}_M \\ \mathbf{G}_M^{\prime} & \mathbf{M} \end{bmatrix}}
#'
#' \strong{Selection response:}
#' \deqn{R_W = k_I \sqrt{\boldsymbol{\beta}_W^{\prime}\mathbf{Q}\boldsymbol{\beta}_W}}
#'
#' \strong{Expected genetic gain per trait:}
#' \deqn{\mathbf{E}_W = k_I \frac{\mathbf{X}\boldsymbol{\beta}_W}{\sqrt{\boldsymbol{\beta}_W^{\prime}\mathbf{Q}\boldsymbol{\beta}_W}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 8.3.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate marker data
#' N_markers <- 100
#' n_traits <- nrow(gmat)
#' G_M <- matrix(rnorm(n_traits * N_markers, sd = 0.5), n_traits, N_markers)
#' M <- diag(N_markers) + matrix(rnorm(N_markers^2, sd = 0.1), N_markers, N_markers)
#' M <- (M + t(M)) / 2 # Make symmetric
#'
#' result <- gw_esim(pmat, gmat, G_M, M)
#' print(result)
#' }
gw_esim <- function(pmat, gmat, G_M, M, selection_intensity = 2.063) {
  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  G_M <- as.matrix(G_M)
  M <- as.matrix(M)

  n_traits <- nrow(pmat)
  n_markers <- ncol(G_M)

  if (!isSymmetric(unname(pmat), tol = 1e-8)) {
    stop("pmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(gmat), tol = 1e-8)) {
    stop("gmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(M), tol = 1e-8)) {
    stop("M must be a symmetric matrix.")
  }
  if (nrow(pmat) != nrow(gmat)) {
    stop("pmat and gmat must have the same dimensions.")
  }
  if (nrow(G_M) != n_traits) {
    stop("G_M must have n_traits rows.")
  }
  if (nrow(M) != n_markers || ncol(M) != n_markers) {
    stop("M must be a square matrix with dimension equal to number of markers.")
  }
  if (n_traits < 2) {
    stop("At least 2 traits are required for GW-ESIM.")
  }

  trait_names <- colnames(pmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  # --------------------------------------------------------------------------
  # Step 1: Construct combined matrices Q and X
  # --------------------------------------------------------------------------
  # Q = [P       G_M ]
  #     [G_M'    M   ]
  Q <- rbind(
    cbind(pmat, G_M),
    cbind(t(G_M), M)
  )

  # X = [C       G_M ]
  #     [G_M'    M   ]
  X <- rbind(
    cbind(gmat, G_M),
    cbind(t(G_M), M)
  )

  # --------------------------------------------------------------------------
  # Step 2: Solve eigenproblem Q^{-1} X
  # --------------------------------------------------------------------------
  Q_inv_X <- .gesim_solve_sym_multi(Q, X)

  ev_result <- .gesim_leading_eigenvector(Q_inv_X)
  lambda2 <- ev_result$value
  b_W <- ev_result$vector

  # --------------------------------------------------------------------------
  # Step 3: Compute metrics
  # --------------------------------------------------------------------------
  metrics <- .genomic_eigen_index_metrics(b_W, Q, X,
    lambda2 = lambda2,
    k_I = selection_intensity
  )

  # Extract corrected eigenvector (sign-corrected for positive quadratic form)
  b_W <- metrics$b

  # Split into phenotype and marker coefficients
  b_y <- b_W[1:n_traits]
  b_m <- b_W[(n_traits + 1):(n_traits + n_markers)]

  # Expected gains are first n_traits elements
  E_W <- metrics$E_vec[1:n_traits]
  names(E_W) <- trait_names

  # --------------------------------------------------------------------------
  # Step 4: Build summary data frame
  # --------------------------------------------------------------------------
  summary_df <- data.frame(
    n_traits = n_traits,
    n_markers = n_markers,
    hI2 = round(metrics$hI2, 6),
    rHI = round(metrics$rHI, 6),
    sigma_I = round(metrics$sigma_I, 6),
    R_W = round(metrics$R, 6),
    lambda2 = round(lambda2, 6),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------------------------
  # Step 5: Return result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    b_y                 = setNames(b_y, trait_names),
    b_m                 = b_m,
    b_combined          = b_W,
    E_W                 = E_W,
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    R_W                 = metrics$R,
    lambda2             = lambda2,
    selection_intensity = selection_intensity,
    Q                   = Q,
    X                   = X,
    trait_names         = trait_names,
    n_markers           = n_markers
  )

  class(result) <- c("gw_esim", "genomic_eigen_index", "list")
  result
}

# ==============================================================================
# 8.4  RGESIM - Restricted Linear Genomic Eigen Selection Index Method
# ==============================================================================

#' Restricted Linear Genomic Eigen Selection Index Method (RGESIM)
#'
#' @description
#' Implements the RGESIM which extends GESIM to allow restrictions on genetic
#' gains of certain traits. Uses the eigen approach with Lagrange multipliers.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param Gamma Covariance between phenotypes and GEBVs (n_traits x n_traits).
#' @param U_mat Restriction matrix (r x n_traits) where r is number of restrictions.
#'   Each row specifies a linear combination of traits to be held at zero gain.
#' @param selection_intensity Selection intensity constant \eqn{k_I}
#'   (default: 2.063 for 10\% selection).
#'
#' @return Object of class \code{"rgesim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with coefficients and metrics.}
#'   \item{\code{b_y}}{Coefficients for phenotypic data.}
#'   \item{\code{b_gamma}}{Coefficients for GEBVs.}
#'   \item{\code{b_combined}}{Combined coefficient vector.}
#'   \item{\code{E_RG}}{Expected genetic gains per trait.}
#'   \item{\code{constrained_response}}{U' * E (should be near zero).}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{hI2}}{Index heritability.}
#'   \item{\code{rHI}}{Accuracy.}
#'   \item{\code{R_RG}}{Selection response.}
#'   \item{\code{lambda2}}{Leading eigenvalue.}
#'   \item{\code{implied_w}}{Implied economic weights.}
#'   \item{\code{K_RG}}{Projection matrix.}
#'   \item{\code{Q_RG}}{Constraint projection matrix.}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Eigenproblem (Section 8.4):}
#' \deqn{(\mathbf{K}_{RG}\mathbf{\Phi}^{-1}\mathbf{A} - \lambda_{RG}^2 \mathbf{I}_{2t})\boldsymbol{\beta}_{RG} = 0}
#'
#' where:
#' \deqn{\mathbf{K}_{RG} = \mathbf{I}_{2t} - \mathbf{Q}_{RG}}
#' \deqn{\mathbf{Q}_{RG} = \mathbf{\Phi}^{-1}\mathbf{A}\mathbf{U}_G(\mathbf{U}_G^{\prime}\mathbf{A}\mathbf{\Phi}^{-1}\mathbf{A}\mathbf{U}_G)^{-1}\mathbf{U}_G^{\prime}\mathbf{A}}
#'
#' \strong{Implied economic weights:}
#' \deqn{\mathbf{w}_{RG} = \mathbf{A}^{-1}[\mathbf{\Phi} + \mathbf{Q}_{RG}^{\prime}\mathbf{A}]\boldsymbol{\beta}_{RG}}
#'
#' \strong{Selection response:}
#' \deqn{R_{RG} = k_I \sqrt{\boldsymbol{\beta}_{RG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}_{RG}}}
#'
#' \strong{Expected genetic gain per trait:}
#' \deqn{\mathbf{E}_{RG} = k_I \frac{\mathbf{A}\boldsymbol{\beta}_{RG}}{\sqrt{\boldsymbol{\beta}_{RG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}_{RG}}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 8.4.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate GEBV covariance
#' Gamma <- gmat * 0.8
#'
#' # Restrict first trait to zero gain
#' # U_mat must be (n_traits x n_restrictions)
#' n_traits <- nrow(gmat)
#' U_mat <- matrix(0, n_traits, 1)
#' U_mat[1, 1] <- 1 # Restrict trait 1
#'
#' result <- rgesim(pmat, gmat, Gamma, U_mat)
#' print(result)
#' print(result$constrained_response) # Should be near zero
#' }
rgesim <- function(pmat, gmat, Gamma, U_mat, selection_intensity = 2.063) {
  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  Gamma <- as.matrix(Gamma)
  U_mat <- as.matrix(U_mat)

  n_traits <- nrow(pmat)
  n_restrictions <- nrow(U_mat)

  if (!isSymmetric(unname(pmat), tol = 1e-8)) {
    stop("pmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(gmat), tol = 1e-8)) {
    stop("gmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(Gamma), tol = 1e-8)) {
    stop("Gamma must be a symmetric matrix.")
  }
  if (nrow(pmat) != nrow(gmat) || nrow(pmat) != nrow(Gamma)) {
    stop("All matrices must have the same dimensions.")
  }
  if (ncol(U_mat) != n_traits) {
    stop("U_mat must have n_traits columns.")
  }
  if (n_traits < 2) {
    stop("At least 2 traits are required for RGESIM.")
  }

  trait_names <- colnames(pmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  # --------------------------------------------------------------------------
  # Step 1: Construct combined matrices Phi and A (same as GESIM)
  # --------------------------------------------------------------------------
  Phi <- rbind(
    cbind(pmat, Gamma),
    cbind(Gamma, Gamma)
  )

  A <- rbind(
    cbind(gmat, Gamma),
    cbind(Gamma, Gamma)
  )

  # --------------------------------------------------------------------------
  # CRITICAL: Implement "2r restriction rule" from Chapter 8.4
  # --------------------------------------------------------------------------
  # When restricting trait gain, we must restrict BOTH:
  #   1. The phenotype coefficient b_y[i] (row i)
  #   2. The GEBV coefficient b_gamma[i] (row t+i)
  # Otherwise, the index bypasses the restriction by shifting weight to GEBV.
  #
  # Solution: Concatenate U_mat horizontally to create (r x 2t) restriction matrix
  # This applies each constraint to both phenotype and GEBV coefficients
  U_G <- cbind(U_mat, U_mat) # [U_mat, U_mat] ensures both y and gamma are restricted

  # --------------------------------------------------------------------------
  # Step 2: Compute constraint projection matrix Q_RG
  # Q_RG = Phi^{-1} A U_G (U_G' A Phi^{-1} A U_G)^{-1} U_G' A
  # --------------------------------------------------------------------------
  Phi_inv_A <- .gesim_solve_sym_multi(Phi, A)

  # Phi^{-1} A U_G
  Phi_inv_A_UG <- Phi_inv_A %*% t(U_G)

  # U_G' A Phi^{-1} A U_G
  middle <- U_G %*% A %*% Phi_inv_A_UG
  middle_inv <- ginv(middle)

  # Complete Q_RG
  Q_RG <- Phi_inv_A_UG %*% middle_inv %*% U_G %*% A


  K_RG <- diag(2 * n_traits) - Q_RG

  # --------------------------------------------------------------------------
  # Step 3: Solve restricted eigenproblem
  # --------------------------------------------------------------------------
  K_RG_Phi_inv_A <- K_RG %*% Phi_inv_A

  ev_result <- .gesim_leading_eigenvector(K_RG_Phi_inv_A)
  lambda2 <- ev_result$value
  b_RG <- ev_result$vector

  # --------------------------------------------------------------------------
  # Step 4: Compute metrics
  # --------------------------------------------------------------------------
  metrics <- .genomic_eigen_index_metrics(b_RG, Phi, A,
    lambda2 = lambda2,
    k_I = selection_intensity
  )

  # Extract corrected eigenvector (sign-corrected for positive quadratic form)
  b_RG <- metrics$b

  # Split coefficients
  b_y <- b_RG[1:n_traits]
  b_gamma <- b_RG[(n_traits + 1):(2 * n_traits)]

  E_RG <- metrics$E_vec[1:n_traits]
  names(E_RG) <- trait_names

  # Verify constraints: U' E should be near zero
  constrained_response <- as.vector(U_mat %*% E_RG)

  # --------------------------------------------------------------------------
  # Step 5: Implied economic weights
  # w_RG = A^{-1} [Phi + Q_RG' A] beta_RG
  # --------------------------------------------------------------------------
  implied_w <- tryCatch(
    {
      A_inv <- ginv(A)
      w_full <- A_inv %*% ((Phi + t(Q_RG) %*% A) %*% b_RG)
      as.numeric(w_full[1:n_traits])
    },
    error = function(e) {
      warning("Could not compute implied economic weights: ", e$message)
      rep(NA_real_, n_traits)
    }
  )
  names(implied_w) <- trait_names

  # --------------------------------------------------------------------------
  # Step 6: Build summary data frame
  # --------------------------------------------------------------------------
  b_y_vec <- round(b_y, 6)
  b_gamma_vec <- round(b_gamma, 6)

  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_len(n_traits))

  b_gamma_df <- as.data.frame(matrix(b_gamma_vec, nrow = 1))
  colnames(b_gamma_df) <- paste0("b_gamma.", seq_len(n_traits))

  summary_df <- data.frame(
    b_y_df,
    b_gamma_df,
    hI2 = round(metrics$hI2, 6),
    rHI = round(metrics$rHI, 6),
    sigma_I = round(metrics$sigma_I, 6),
    R_RG = round(metrics$R, 6),
    lambda2 = round(lambda2, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # --------------------------------------------------------------------------
  # Step 7: Return result
  # --------------------------------------------------------------------------
  result <- list(
    summary               = summary_df,
    b_y                   = setNames(b_y, trait_names),
    b_gamma               = setNames(b_gamma, trait_names),
    b_combined            = b_RG,
    E_RG                  = E_RG,
    constrained_response  = constrained_response,
    sigma_I               = metrics$sigma_I,
    hI2                   = metrics$hI2,
    rHI                   = metrics$rHI,
    R_RG                  = metrics$R,
    lambda2               = lambda2,
    implied_w             = implied_w,
    K_RG                  = K_RG,
    Q_RG                  = Q_RG,
    U_mat                 = U_mat,
    selection_intensity   = selection_intensity,
    Phi                   = Phi,
    A                     = A,
    trait_names           = trait_names,
    n_restrictions        = n_restrictions
  )

  class(result) <- c("rgesim", "genomic_eigen_index", "list")
  result
}

# ==============================================================================
# 8.5  PPG-GESIM - Predetermined Proportional Gain Genomic Eigen Selection Index
# ==============================================================================

#' Predetermined Proportional Gain Genomic Eigen Selection Index (PPG-GESIM)
#'
#' @description
#' Implements the PPG-GESIM which extends GESIM to enforce that genetic gains
#' are proportional to a user-specified vector d. Combines eigen approach with
#' predetermined gain proportions.
#'
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param Gamma Covariance between phenotypes and GEBVs (n_traits x n_traits).
#' @param d Numeric vector of desired proportional gains (length n_traits).
#'   The ratios among elements define target gain proportions.
#' @param selection_intensity Selection intensity constant \eqn{k_I}
#'   (default: 2.063 for 10\% selection).
#'
#' @return Object of class \code{"ppg_gesim"}, a list with:
#' \describe{
#'   \item{\code{summary}}{Data frame with coefficients and metrics.}
#'   \item{\code{b_y}}{Coefficients for phenotypic data.}
#'   \item{\code{b_gamma}}{Coefficients for GEBVs.}
#'   \item{\code{b_combined}}{Combined coefficient vector.}
#'   \item{\code{E_PG}}{Expected genetic gains per trait.}
#'   \item{\code{gain_ratios}}{Ratios of actual to desired gains (should be constant).}
#'   \item{\code{d}}{Original desired proportional gains (length t).}
#'   \item{\code{d_PG}}{Extended proportional gains (length 2t, includes GEBV targets).}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{hI2}}{Index heritability.}
#'   \item{\code{rHI}}{Accuracy.}
#'   \item{\code{R_PG}}{Selection response.}
#'   \item{\code{lambda2}}{Leading eigenvalue.}
#'   \item{\code{implied_w}}{Implied economic weights.}
#'   \item{\code{U_PG}}{Restriction matrix ((2t-1) x 2t).}
#'   \item{\code{selection_intensity}}{Selection intensity used.}
#' }
#'
#' @details
#' \strong{Eigenproblem (Section 8.5):}
#' \deqn{(\mathbf{T}_{PG} - \lambda_{PG}^2 \mathbf{I}_{2t})\boldsymbol{\beta}_{PG} = 0}
#'
#' where:
#' \deqn{\mathbf{T}_{PG} = \mathbf{K}_{RG}\mathbf{\Phi}^{-1}\mathbf{A} + \mathbf{B}}
#' \deqn{\mathbf{B} = \boldsymbol{\delta}\boldsymbol{\varphi}^{\prime}}
#'
#' \strong{Implied economic weights:}
#' \deqn{\mathbf{w}_{PG} = \mathbf{A}^{-1}[\mathbf{\Phi} + \mathbf{Q}_{PG}^{\prime}\mathbf{A}]\boldsymbol{\beta}_{PG}}
#'
#' \strong{Selection response:}
#' \deqn{R_{PG} = k_I \sqrt{\boldsymbol{\beta}_{PG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}_{PG}}}
#'
#' \strong{Expected genetic gain per trait:}
#' \deqn{\mathbf{E}_{PG} = k_I \frac{\mathbf{A}\boldsymbol{\beta}_{PG}}{\sqrt{\boldsymbol{\beta}_{PG}^{\prime}\mathbf{\Phi}\boldsymbol{\beta}_{PG}}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Section 8.5.
#'
#' @export
#' @examples
#' \dontrun{
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate GEBV covariance
#' Gamma <- gmat * 0.8
#'
#' # Desired proportional gains (e.g., 2:1:3 ratio for first 3 traits)
#' d <- c(2, 1, 3, 1, 1, 1, 1)
#'
#' result <- ppg_gesim(pmat, gmat, Gamma, d)
#' print(result)
#' print(result$gain_ratios) # Should be approximately constant
#' }
ppg_gesim <- function(pmat, gmat, Gamma, d, selection_intensity = 2.063) {
  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)
  Gamma <- as.matrix(Gamma)
  d <- as.numeric(d)

  n_traits <- nrow(pmat)

  if (!isSymmetric(unname(pmat), tol = 1e-8)) {
    stop("pmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(gmat), tol = 1e-8)) {
    stop("gmat must be a symmetric matrix.")
  }
  if (!isSymmetric(unname(Gamma), tol = 1e-8)) {
    stop("Gamma must be a symmetric matrix.")
  }
  if (nrow(pmat) != nrow(gmat) || nrow(pmat) != nrow(Gamma)) {
    stop("All matrices must have the same dimensions.")
  }
  if (length(d) != n_traits) {
    stop("d must have length equal to number of traits.")
  }
  if (n_traits < 2) {
    stop("At least 2 traits are required for PPG-GESIM.")
  }

  trait_names <- colnames(pmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait_", seq_len(n_traits))
  }

  # --------------------------------------------------------------------------
  # CRITICAL: Extend d to d_PG for genomic space (Chapter 8.5)
  # --------------------------------------------------------------------------
  # Section 8.5 states: "the vector of PPG (d_PG) should be twice the
  # standard vector... d_PG' = [d_1...d_t, d_t+1...d_2t]"
  #
  # Since Phi and A are 2tx2t (phenotypes + GEBVs), we need a 2t-length
  # proportional gain vector. Assuming equal targets for phenotype and GEBV:
  d_PG <- c(d, d) # [phenotype targets, GEBV targets]

  # Note: d_PG now has length 2t, matching the dimension of beta = [beta_y; beta_gamma]

  # --------------------------------------------------------------------------
  # Step 1: Construct combined matrices Phi and A
  # --------------------------------------------------------------------------
  Phi <- rbind(
    cbind(pmat, Gamma),
    cbind(Gamma, Gamma)
  )

  A <- rbind(
    cbind(gmat, Gamma),
    cbind(Gamma, Gamma)
  )

  # --------------------------------------------------------------------------
  # Step 2: Construct restriction matrix from d_PG (now 2t-length)
  # Build (2t-1) restrictions: d_PG[i+1] * g[i] - d_PG[i] * g[i+1] = 0
  # --------------------------------------------------------------------------
  n_combined <- 2 * n_traits # Work in 2t space
  U_PG <- matrix(0, n_combined - 1, n_combined)
  for (i in seq_len(n_combined - 1)) {
    U_PG[i, i] <- d_PG[i + 1]
    U_PG[i, i + 1] <- -d_PG[i]
  }

  # --------------------------------------------------------------------------
  # Step 3: Compute projection matrices (same structure as RGESIM)
  # --------------------------------------------------------------------------
  # Now U_PG is (2t-1) x 2t, matching the dimension of Phi and A (2t x 2t)
  Phi_inv_A <- .gesim_solve_sym_multi(Phi, A)

  Phi_inv_A_UPG <- Phi_inv_A %*% t(U_PG)
  middle <- U_PG %*% A %*% Phi_inv_A_UPG
  middle_inv <- ginv(middle)

  Q_PG <- Phi_inv_A_UPG %*% middle_inv %*% U_PG %*% A
  K_PG <- diag(2 * n_traits) - Q_PG

  # --------------------------------------------------------------------------
  # Step 4: Construct matrix B = delta %*% t(phi)
  # delta = Phi^{-1} A d_PG  (now using the full 2t-length vector)
  # phi is d_PG (the 2t-length proportional gain vector)
  # --------------------------------------------------------------------------
  delta <- Phi_inv_A %*% d_PG
  phi <- d_PG


  B <- delta %*% t(phi)


  T_PG <- K_PG %*% Phi_inv_A + B

  # --------------------------------------------------------------------------
  # Step 5: Solve eigenproblem for T_PG
  # --------------------------------------------------------------------------
  ev_result <- .gesim_leading_eigenvector(T_PG)
  lambda2 <- ev_result$value
  b_PG <- ev_result$vector

  # --------------------------------------------------------------------------
  # Step 6: Compute metrics
  # --------------------------------------------------------------------------
  metrics <- .genomic_eigen_index_metrics(b_PG, Phi, A,
    lambda2 = lambda2,
    k_I = selection_intensity
  )

  # Extract corrected eigenvector (sign-corrected for positive quadratic form)
  b_PG <- metrics$b

  # Split coefficients
  b_y <- b_PG[1:n_traits]
  b_gamma <- b_PG[(n_traits + 1):(2 * n_traits)]

  E_PG <- metrics$E_vec[1:n_traits]
  names(E_PG) <- trait_names

  # Compute gain ratios: E_PG / d (should be approximately constant)
  gain_ratios <- E_PG / d
  gain_ratios[!is.finite(gain_ratios)] <- NA_real_

  # --------------------------------------------------------------------------
  # Step 7: Implied economic weights
  # --------------------------------------------------------------------------
  implied_w <- tryCatch(
    {
      A_inv <- ginv(A)
      w_full <- A_inv %*% ((Phi + t(Q_PG) %*% A) %*% b_PG)
      as.numeric(w_full[1:n_traits])
    },
    error = function(e) {
      warning("Could not compute implied economic weights: ", e$message)
      rep(NA_real_, n_traits)
    }
  )
  names(implied_w) <- trait_names

  # --------------------------------------------------------------------------
  # Step 8: Build summary data frame
  # --------------------------------------------------------------------------
  b_y_vec <- round(b_y, 6)
  b_gamma_vec <- round(b_gamma, 6)

  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_len(n_traits))

  b_gamma_df <- as.data.frame(matrix(b_gamma_vec, nrow = 1))
  colnames(b_gamma_df) <- paste0("b_gamma.", seq_len(n_traits))

  summary_df <- data.frame(
    b_y_df,
    b_gamma_df,
    hI2 = round(metrics$hI2, 6),
    rHI = round(metrics$rHI, 6),
    sigma_I = round(metrics$sigma_I, 6),
    R_PG = round(metrics$R, 6),
    lambda2 = round(lambda2, 6),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # --------------------------------------------------------------------------
  # Step 9: Return result
  # --------------------------------------------------------------------------
  result <- list(
    summary             = summary_df,
    b_y                 = setNames(b_y, trait_names),
    b_gamma             = setNames(b_gamma, trait_names),
    b_combined          = b_PG,
    E_PG                = E_PG,
    gain_ratios         = gain_ratios,
    d                   = d,
    d_PG                = d_PG, # Store the expanded 2t-length vector
    sigma_I             = metrics$sigma_I,
    hI2                 = metrics$hI2,
    rHI                 = metrics$rHI,
    R_PG                = metrics$R,
    lambda2             = lambda2,
    implied_w           = implied_w,
    K_PG                = K_PG,
    Q_PG                = Q_PG,
    U_PG                = U_PG, # Store restriction matrix (was U_mat)
    selection_intensity = selection_intensity,
    Phi                 = Phi,
    A                   = A,
    trait_names         = trait_names
  )

  class(result) <- c("ppg_gesim", "genomic_eigen_index", "list")
  result
}

# ==============================================================================
# S3 PRINT METHODS FOR GENOMIC EIGEN INDICES
# ==============================================================================

#' Print method for MESIM
#' @param x Object of class 'mesim'
#' @param ... Additional arguments (unused)
#' @export
print.mesim <- function(x, ...) {
  cat("\n==============================================================\n")
  cat("MOLECULAR EIGEN SELECTION INDEX METHOD (MESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.1\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_M^2 (h^2_I):     %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy (r_HI):        %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev (sigma_I): %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response (R_M): %.6f\n", x$R_M))

  cat("\n-------------------------------------------------------------\n")
  cat("PHENOTYPE COEFFICIENTS (b_y)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_y = round(x$b_y, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("MARKER SCORE COEFFICIENTS (b_s)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_s = round(x$b_s, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_M)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, E_M = round(x$E_M, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' Print method for GESIM
#' @param x Object of class 'gesim'
#' @param ... Additional arguments (unused)
#' @export
print.gesim <- function(x, ...) {
  cat("\n==============================================================\n")
  cat("LINEAR GENOMIC EIGEN SELECTION INDEX METHOD (GESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.2\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_G^2 (h^2_I):     %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy (r_HI):        %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev (sigma_I): %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response (R_G): %.6f\n", x$R_G))

  cat("\n-------------------------------------------------------------\n")
  cat("PHENOTYPE COEFFICIENTS (b_y)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_y = round(x$b_y, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("GEBV COEFFICIENTS (b_gamma)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_gamma = round(x$b_gamma, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_G)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, E_G = round(x$E_G, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("IMPLIED ECONOMIC WEIGHTS (w_G)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, Implied_w = round(x$implied_w, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' Print method for GW-ESIM
#' @param x Object of class 'gw_esim'
#' @param ... Additional arguments (unused)
#' @export
print.gw_esim <- function(x, ...) {
  cat("\n==============================================================\n")
  cat("GENOME-WIDE LINEAR EIGEN SELECTION INDEX METHOD (GW-ESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.3\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n")
  cat("Number of markers:        ", x$n_markers, "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_W^2 (h^2_I):     %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy (r_HI):        %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev (sigma_I): %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response (R_W): %.6f\n", x$R_W))

  cat("\n-------------------------------------------------------------\n")
  cat("PHENOTYPE COEFFICIENTS (b_y)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_y = round(x$b_y, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("MARKER COEFFICIENTS (b_m) - First 10\n")
  cat("-------------------------------------------------------------\n")
  n_show <- min(10, length(x$b_m))
  print(data.frame(
    Marker = 1:n_show, b_m = round(x$b_m[1:n_show], 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)
  if (length(x$b_m) > 10) {
    cat("  ... and", length(x$b_m) - 10, "more markers\n")
  }

  cat("\n-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_W)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, E_W = round(x$E_W, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' Print method for RGESIM
#' @param x Object of class 'rgesim'
#' @param ... Additional arguments (unused)
#' @export
print.rgesim <- function(x, ...) {
  cat("\n==============================================================\n")
  cat("RESTRICTED LINEAR GENOMIC EIGEN SELECTION INDEX (RGESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.4\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n")
  cat("Number of restrictions:   ", x$n_restrictions, "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_RG^2 (h^2_I):    %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy (r_HI):        %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev (sigma_I): %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response (R_RG): %.6f\n", x$R_RG))

  cat("\n-------------------------------------------------------------\n")
  cat("PHENOTYPE COEFFICIENTS (b_y)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_y = round(x$b_y, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("GEBV COEFFICIENTS (b_gamma)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_gamma = round(x$b_gamma, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_RG)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, E_RG = round(x$E_RG, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("CONSTRAINT VERIFICATION\n")
  cat("-------------------------------------------------------------\n")
  cat("Constrained response (should be near zero):\n")
  print(round(x$constrained_response, 8))

  cat("\n-------------------------------------------------------------\n")
  cat("IMPLIED ECONOMIC WEIGHTS (w_RG)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, Implied_w = round(x$implied_w, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' Print method for PPG-GESIM
#' @param x Object of class 'ppg_gesim'
#' @param ... Additional arguments (unused)
#' @export
print.ppg_gesim <- function(x, ...) {
  cat("\n==============================================================\n")
  cat("PREDETERMINED PROPORTIONAL GAIN GENOMIC EIGEN INDEX (PPG-GESIM)\n")
  cat("Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.5\n")
  cat("==============================================================\n\n")

  cat("Selection intensity (k_I):", x$selection_intensity, "\n")
  cat("Number of traits:         ", length(x$trait_names), "\n\n")

  cat("-------------------------------------------------------------\n")
  cat("DESIRED PROPORTIONAL GAINS (d)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, d = round(x$d, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf("  lambda_PG^2 (h^2_I):    %.6f\n", x$lambda2))
  cat(sprintf("  Accuracy (r_HI):        %.6f\n", x$rHI))
  cat(sprintf("  Index Std Dev (sigma_I): %.6f\n", x$sigma_I))
  cat(sprintf("  Selection Response (R_PG): %.6f\n", x$R_PG))

  cat("\n-------------------------------------------------------------\n")
  cat("PHENOTYPE COEFFICIENTS (b_y)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_y = round(x$b_y, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("GEBV COEFFICIENTS (b_gamma)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, b_gamma = round(x$b_gamma, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT (E_PG)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, E_PG = round(x$E_PG, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n-------------------------------------------------------------\n")
  cat("PROPORTIONAL GAIN VERIFICATION\n")
  cat("-------------------------------------------------------------\n")
  cat("Gain ratios (E_PG / d) - should be approximately constant:\n")
  print(data.frame(
    Trait = x$trait_names, Ratio = round(x$gain_ratios, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)
  cat("Standard deviation of ratios:", round(sd(x$gain_ratios, na.rm = TRUE), 6), "\n")

  cat("\n-------------------------------------------------------------\n")
  cat("IMPLIED ECONOMIC WEIGHTS (w_PG)\n")
  cat("-------------------------------------------------------------\n")
  print(data.frame(
    Trait = x$trait_names, Implied_w = round(x$implied_w, 6),
    stringsAsFactors = FALSE
  ), row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' Summary method for MESIM
#' @param object Object of class 'mesim'
#' @param ... Additional arguments (unused)
#' @export
summary.mesim <- function(object, ...) {
  print(object, ...)
  cat("==============================================================\n")
  cat("SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' Summary method for GESIM
#' @param object Object of class 'gesim'
#' @param ... Additional arguments (unused)
#' @export
summary.gesim <- function(object, ...) {
  print(object, ...)
  cat("==============================================================\n")
  cat("SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' Summary method for GW-ESIM
#' @param object Object of class 'gw_esim'
#' @param ... Additional arguments (unused)
#' @export
summary.gw_esim <- function(object, ...) {
  print(object, ...)
  cat("==============================================================\n")
  cat("SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' Summary method for RGESIM
#' @param object Object of class 'rgesim'
#' @param ... Additional arguments (unused)
#' @export
summary.rgesim <- function(object, ...) {
  print(object, ...)
  cat("==============================================================\n")
  cat("SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")
  invisible(object)
}

#' Summary method for PPG-GESIM
#' @param object Object of class 'ppg_gesim'
#' @param ... Additional arguments (unused)
#' @export
summary.ppg_gesim <- function(object, ...) {
  print(object, ...)
  cat("==============================================================\n")
  cat("SUMMARY TABLE\n")
  cat("==============================================================\n\n")
  print(object$summary, row.names = FALSE)
  cat("\n")
  invisible(object)
}
