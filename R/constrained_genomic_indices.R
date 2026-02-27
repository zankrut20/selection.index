#' Constrained Genomic Selection Indices
#' @name constrained_genomic_indices
#'
#' @description
#' Implements constrained genomic selection index methods using GEBVs.
#' These methods allow breeders to impose restrictions on genetic gains
#' while maintaining genomic selection efficiency.
#'
#' Methods included:
#' - Restricted Linear Genomic Selection Index (RLGSI)
#' - Predetermined Proportional Gains Linear Genomic Selection Index (PPG-LGSI)
#' - Combined Restricted Linear Genomic Selection Index (CRLGSI)
#' - Combined Predetermined Proportional Gains Linear Genomic Selection Index (CPPG-LGSI)
#'
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 6.
#'
#' @keywords internal
#' @importFrom stats cov
#' @importFrom MASS ginv
NULL

# ==============================================================================
# HELPER FUNCTIONS (using C++ primitives for all math operations)
# ==============================================================================

#' Compute genomic index metrics using C++ primitives
#' @keywords internal
#' @noRd
.genomic_index_metrics <- function(b, Gamma, w = NULL, k_I = 2.063, L_G = 1, GAY = NULL, gmat = NULL) {
  b <- as.numeric(b)

  # Use C++ primitives for all calculations
  bGb <- cpp_quadratic_form_sym(b, Gamma)
  sigma_I <- if (bGb > 0) sqrt(bGb) else NA_real_

  # Selection response: R = (k_I / L_G) * sqrt(b'Γb)
  R <- if (!is.na(sigma_I)) (k_I / L_G) * sigma_I else NA_real_

  # Expected genetic gain per trait: E = (k_I / L_G) * (Γb) / sqrt(b'Γb)
  E_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    (k_I / L_G) * as.vector(Gamma %*% b) / sigma_I
  } else {
    rep(NA_real_, nrow(Gamma))
  }

  # Overall genetic advance (if weights provided)
  GA <- NA_real_
  PRE <- NA_real_
  if (!is.null(w)) {
    bGw <- cpp_quadratic_form(b, Gamma, w)
    GA <- if (!is.na(sigma_I) && sigma_I > 0) (k_I / L_G) * bGw / sigma_I else NA_real_
    PRE_constant <- if (is.null(GAY)) 100 else 100 / GAY
    PRE <- if (!is.na(GA)) GA * PRE_constant else NA_real_
  }

  # Index accuracy (for genomic indices, this is the correlation with breeding objective)
  rHI <- if (!is.null(w)) {
    # Use true genetic variance if provided, otherwise use GEBV variance as approximation
    # Var(H) = w'Cw where C is true genetic variance (gmat)
    C <- if (!is.null(gmat)) gmat else Gamma
    wCw <- cpp_quadratic_form_sym(w, C)
    if (!is.na(bGb) && !is.na(wCw) && bGb > 0 && wCw > 0) {
      bGw_val <- cpp_quadratic_form(b, Gamma, w)
      pmin(pmax(bGw_val / (sqrt(bGb) * sqrt(wCw)), 0), 1)
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }

  list(
    bGb = bGb,
    sigma_I = sigma_I,
    R = R,
    E_vec = as.vector(E_vec),
    GA = GA,
    PRE = PRE,
    rHI = rHI
  )
}

#' Compute combined index metrics (for CRLGSI/CPPG-LGSI)
#' @keywords internal
#' @noRd
.combined_index_metrics <- function(b, T_C, Psi_C, w = NULL, k_I = 2.063, L_I = 1, GAY = NULL) {
  b <- as.numeric(b)

  # Use C++ primitives for all calculations
  bTb <- cpp_quadratic_form_sym(b, T_C)
  sigma_I <- if (bTb > 0) sqrt(bTb) else NA_real_

  # Selection response: R = (k_I / L_I) * sqrt(b'T_Cb)
  R <- if (!is.na(sigma_I)) (k_I / L_I) * sigma_I else NA_real_

  # Expected genetic gain per trait: E = (k_I / L_I) * (Ψ_C'b)[1:t] / sqrt(b'T_Cb)
  # Psi_C is 2t x 2t, b is 2t x 1, Psi_C'b is 2t x 1
  # We extract first t elements (genetic gains for traits)
  n_traits <- ncol(Psi_C) / 2
  E_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    Psi_b <- as.vector(crossprod(Psi_C, b))
    (k_I / L_I) * Psi_b[1:n_traits] / sigma_I
  } else {
    rep(NA_real_, n_traits)
  }

  # Overall genetic advance (if weights provided)
  GA <- NA_real_
  PRE <- NA_real_
  if (!is.null(w)) {
    # a_C = [w; 0] - weights apply only to breeding values, not GEBVs
    n_traits <- length(w)
    a_C <- c(w, rep(0, n_traits))
    bPsiw <- as.numeric(crossprod(b, Psi_C) %*% a_C)
    GA <- if (!is.na(sigma_I) && sigma_I > 0) (k_I / L_I) * bPsiw / sigma_I else NA_real_
    PRE_constant <- if (is.null(GAY)) 100 else 100 / GAY
    PRE <- if (!is.na(GA)) GA * PRE_constant else NA_real_
  }

  # Index accuracy
  rHI <- if (!is.null(w)) {
    # For combined index, H = w'g, Var(H) = w'Gw where G = Psi_C[1:t, 1:t]
    n_traits <- length(w)
    G <- Psi_C[1:n_traits, 1:n_traits, drop = FALSE]
    wGw <- cpp_quadratic_form_sym(w, G)
    if (!is.na(bTb) && !is.na(wGw) && bTb > 0 && wGw > 0) {
      # a_C = [w; 0] - weights apply only to breeding values
      a_C <- c(w, rep(0, n_traits))
      bPsiw_val <- as.numeric(crossprod(b, Psi_C) %*% a_C)
      pmin(pmax(bPsiw_val / (sqrt(bTb) * sqrt(wGw)), 0), 1)
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }

  list(
    bTb = bTb,
    sigma_I = sigma_I,
    R = R,
    E_vec = as.vector(E_vec),
    GA = GA,
    PRE = PRE,
    rHI = rHI
  )
}

# ==============================================================================
# 6.1 RESTRICTED LINEAR GENOMIC SELECTION INDEX (RLGSI)
# ==============================================================================

#' Restricted Linear Genomic Selection Index (RLGSI)
#'
#' @description
#' Implements the Restricted Linear Genomic Selection Index where genetic gains
#' are constrained to zero for specific traits while maximizing gains for others.
#' Uses GEBVs only (no phenotypic data required).
#'
#' @param Gamma GEBV variance-covariance matrix (n_traits x n_traits).
#'   This represents the variance of GEBVs, typically computed from predicted breeding values.
#' @param wmat Economic weights matrix (n_traits x k), or vector
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param restricted_traits Vector of trait indices to restrict (default: NULL).
#'   Example: c(1, 3) restricts traits 1 and 3 to zero gain.
#' @param U Constraint matrix (n_traits x n_constraints). Each column defines a restriction.
#'   Alternative to restricted_traits for custom constraints. Ignored if restricted_traits is provided.
#' @param k_I Selection intensity (default: 2.063 for 10 percent selection)
#' @param L_G Standardization constant (default: 1). Can be set to sqrt(w'Gw) for standardization.
#' @param gmat Optional. True genetic variance-covariance matrix for exact accuracy calculation.
#'   If NULL, uses Gamma as approximation. Providing gmat ensures textbook-perfect accuracy metric.
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients, response metrics
#'     \item \code{b} - Vector of RLGSI coefficients (β_RG)
#'     \item \code{E} - Named vector of expected genetic gains per trait
#'     \item \code{R} - Overall selection response
#'     \item \code{U} - Constraint matrix used
#'     \item \code{constrained_response} - Realized gains for constrained traits (should be ~0)
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 6, Section 6.1):}
#'
#' The RLGSI minimizes the mean squared difference between the index I = β'γ and
#' the breeding objective H = w'g under the restriction: U'Γβ = 0
#'
#' Solution involves solving the augmented system:
#' \deqn{\begin{bmatrix} \Gamma & \Gamma U \\ U'\Gamma & 0 \end{bmatrix}
#'       \begin{bmatrix} \beta \\ v \end{bmatrix} =
#'       \begin{bmatrix} \Gamma w \\ 0 \end{bmatrix}}
#'
#' Where:
#' - Γ (Gamma) = Var(GEBVs) - GEBV variance-covariance matrix
#' - U = Constraint matrix (each column is a restriction vector)
#' - w = Economic weights
#' - β_RG = RLGSI coefficient vector
#' - v = Lagrange multipliers
#'
#' Selection response: \eqn{R_{RG} = (k_I / L_G) * sqrt(beta_RG' * Gamma * beta_RG)}
#'
#' Expected gains: \eqn{E_{RG} = (k_I / L_G) * (Gamma * beta_RG) / sqrt(beta_RG' * Gamma * beta_RG)}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate GEBV variance-covariance matrix
#' set.seed(123)
#' n_traits <- 5
#' Gamma <- matrix(rnorm(n_traits^2), n_traits, n_traits)
#' Gamma <- (Gamma + t(Gamma)) / 2 # Make symmetric
#' diag(Gamma) <- abs(diag(Gamma)) + 2 # Ensure positive definite
#'
#' # Economic weights
#' w <- c(10, 8, 6, 4, 2)
#'
#' # Restrict traits 2 and 4 to zero gain
#' result <- rlgsi(Gamma, w, restricted_traits = c(2, 4))
#' print(result$summary)
#' print(result$E) # Check that traits 2 and 4 have ~0 gain
#' }
rlgsi <- function(Gamma, wmat, wcol = 1,
                  restricted_traits = NULL, U = NULL,
                  k_I = 2.063, L_G = 1, gmat = NULL, GAY = NULL) {
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  Gamma <- as.matrix(Gamma)
  wmat <- as.matrix(wmat)
  n_traits <- nrow(Gamma)

  if (ncol(Gamma) != n_traits) {
    stop("Gamma must be a square matrix")
  }

  if (nrow(wmat) != n_traits) {
    stop("wmat must have ", n_traits, " rows")
  }

  # Auto-create U from restricted_traits
  if (!is.null(restricted_traits)) {
    if (!is.numeric(restricted_traits) || any(restricted_traits < 1) || any(restricted_traits > n_traits)) {
      stop("restricted_traits must be valid trait indices (1 to ", n_traits, ")")
    }
    U <- diag(n_traits)[, restricted_traits, drop = FALSE]
  } else if (is.null(U)) {
    stop("Either 'restricted_traits' or 'U' must be provided")
  }

  U <- as.matrix(U)
  if (nrow(U) != n_traits) {
    stop("U must have ", n_traits, " rows")
  }

  # Extract weight vector using C++ primitive
  w <- cpp_extract_vector(wmat, seq_len(n_traits), wcol - 1L)

  n_constraints <- ncol(U)

  # ============================================================================
  # SOLVE AUGMENTED SYSTEM (Using R logic, C++ for matrix ops)
  # ============================================================================

  # Build augmented matrix: [Γ, ΓU; U'Γ, 0]
  GammaU <- Gamma %*% U
  UtGamma <- t(U) %*% Gamma
  zeros <- matrix(0, n_constraints, n_constraints)

  augmented_mat <- rbind(
    cbind(Gamma, GammaU),
    cbind(UtGamma, zeros)
  )

  # Build RHS: [Γw; 0]
  Gammaw <- Gamma %*% w
  rhs <- c(Gammaw, rep(0, n_constraints))

  # Solve system using MASS::ginv for robustness
  solution <- ginv(augmented_mat) %*% rhs

  # Extract β_RG (first n_traits elements)
  b_RG <- solution[1:n_traits]

  # Check for numerical issues
  if (any(is.na(b_RG)) || any(is.infinite(b_RG))) {
    stop("RLGSI coefficients contain NA or Inf. Check matrix conditioning.")
  }

  # ============================================================================
  # COMPUTE METRICS
  # ============================================================================

  metrics <- .genomic_index_metrics(b_RG, Gamma, w, k_I, L_G,
    GAY = if (missing(GAY)) NULL else GAY,
    gmat = gmat
  )

  # Check constraint satisfaction: U'Γβ should be ~0
  constrained_response <- as.vector(t(U) %*% Gamma %*% b_RG)

  # ============================================================================
  # BUILD SUMMARY OUTPUT
  # ============================================================================

  b_vec <- as.numeric(b_RG)
  b_vec <- round(b_vec, 4)

  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_along(b_vec)) # seq_len(length(b_vec)))

  summary_df <- data.frame(
    b_df,
    R = round(metrics$R, 4),
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    rHI = round(metrics$rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Name vectors
  trait_names <- colnames(Gamma)
  if (!is.null(trait_names)) {
    names(metrics$E_vec) <- trait_names
  }

  rownames(summary_df) <- NULL

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  result <- list(
    summary = summary_df,
    b = b_vec,
    E = metrics$E_vec,
    R = metrics$R,
    GA = metrics$GA,
    PRE = metrics$PRE,
    rHI = metrics$rHI,
    sigma_I = metrics$sigma_I,
    U = U,
    constrained_response = constrained_response,
    k_I = k_I,
    L_G = L_G
  )

  class(result) <- c("rlgsi", "genomic_index", "list")

  result
}

# ==============================================================================
# 6.2 PREDETERMINED PROPORTIONAL GAINS LINEAR GENOMIC SELECTION INDEX (PPG-LGSI)
# ==============================================================================

#' Predetermined Proportional Gains Linear Genomic Selection Index (PPG-LGSI)
#'
#' @description
#' Implements the PPG-LGSI where breeders specify desired proportional gains
#' between traits rather than restricting specific traits to zero.
#' This is genomic version of PPG-LPSI using GEBVs only.
#'
#' @param Gamma GEBV variance-covariance matrix (n_traits x n_traits)
#' @param d Vector of desired proportional gains (length n_traits or n_constraints).
#'   If length n_traits, constraints are applied to all traits.
#'   If length n_constraints, must provide U matrix.
#' @param wmat Optional. Economic weights for GA/PRE calculation
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param U Optional. Constraint matrix (n_traits x n_constraints).
#'   If NULL, assumes d applies to all traits (U = I).
#' @param k_I Selection intensity (default: 2.063)
#' @param L_G Standardization constant (default: 1)
#' @param gmat Optional. True genetic variance-covariance matrix for exact accuracy calculation.
#'   If NULL, uses Gamma as approximation.
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients and metrics
#'     \item \code{b} - Vector of PPG-LGSI coefficients (β_PG)
#'     \item \code{E} - Named vector of expected genetic gains per trait
#'     \item \code{theta_G} - Proportionality constant
#'     \item \code{gain_ratios} - Ratios of achieved to desired gains
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 6, Section 6.2):}
#'
#' Alternative form: \eqn{beta_PG = beta_RG + theta_G * U * (U' * Gamma * U)^{-1} * d}
#'
#' Where:
#' - beta_RG = Restricted index coefficients (from RLGSI)
#' - theta_G = Proportionality constant
#' - d = Vector of desired proportional gains
#'
#' Proportionality constant:
#' \deqn{theta_G = (d' * (U' * Gamma * U)^{-1} * U' * Gamma * w) / (d' * (U' * Gamma * U)^{-1} * d)}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate GEBV variance-covariance matrix
#' set.seed(123)
#' n_traits <- 5
#' Gamma <- matrix(rnorm(n_traits^2), n_traits, n_traits)
#' Gamma <- (Gamma + t(Gamma)) / 2
#' diag(Gamma) <- abs(diag(Gamma)) + 2
#'
#' # Desired proportional gains (e.g., 2:1:1:0:0 ratio)
#' d <- c(2, 1, 1, 0, 0)
#'
#' # Economic weights
#' w <- c(10, 8, 6, 4, 2)
#'
#' result <- ppg_lgsi(Gamma, d, wmat = w)
#' print(result$summary)
#' print(result$gain_ratios) # Should be approximately proportional to d
#' }
ppg_lgsi <- function(Gamma, d, wmat = NULL, wcol = 1, U = NULL,
                     k_I = 2.063, L_G = 1, gmat = NULL, GAY = NULL) {
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  Gamma <- as.matrix(Gamma)
  d <- as.numeric(d)
  n_traits <- nrow(Gamma)

  if (ncol(Gamma) != n_traits) {
    stop("Gamma must be a square matrix")
  }

  # Handle U matrix
  if (is.null(U)) {
    # If U not provided, assume d applies to all traits
    if (length(d) != n_traits) {
      stop("d must have length ", n_traits, " when U is not provided")
    }
    U <- diag(n_traits)
  } else {
    U <- as.matrix(U)
    if (nrow(U) != n_traits) {
      stop("U must have ", n_traits, " rows")
    }
    if (length(d) != ncol(U)) {
      stop("d must have length ", ncol(U), " (number of constraints)")
    }
  }

  # Handle weights
  w <- NULL
  if (!is.null(wmat)) {
    wmat <- as.matrix(wmat)
    if (nrow(wmat) != n_traits) {
      stop("wmat must have ", n_traits, " rows")
    }
    w <- cpp_extract_vector(wmat, seq_len(n_traits), wcol - 1L)
  }

  # ============================================================================
  # COMPUTE β_RG (RLGSI with zero restrictions)
  # ============================================================================

  # For PPG-LGSI, first compute RLGSI with same U
  # Create weights if not provided (use equal weights for RLGSI component)
  if (is.null(w)) {
    w_temp <- rep(1, n_traits)
  } else {
    w_temp <- w
  }

  # Compute RLGSI component
  rlgsi_result <- rlgsi(Gamma, w_temp,
    wcol = 1, U = U,
    k_I = k_I, L_G = L_G, gmat = gmat, GAY = NULL
  )
  b_RG <- rlgsi_result$b

  # ============================================================================
  # COMPUTE PROPORTIONALITY CONSTANT θ_G
  # ============================================================================

  # U'ΓU
  UtGammaU <- t(U) %*% Gamma %*% U

  # (U'ΓU)^(-1)
  UtGammaU_inv <- solve(UtGammaU)

  # θ_G numerator: d'(U'ΓU)^(-1)U'Γw
  if (!is.null(w)) {
    UtGammaw <- t(U) %*% Gamma %*% w
    theta_num <- as.numeric(t(d) %*% UtGammaU_inv %*% UtGammaw)
  } else {
    theta_num <- 0
  }

  # θ_G denominator: d'(U'ΓU)^(-1)d
  theta_denom <- as.numeric(t(d) %*% UtGammaU_inv %*% d)

  theta_G <- if (abs(theta_denom) > 1e-10) theta_num / theta_denom else 0

  # ============================================================================
  # COMPUTE β_PG = β_RG + θ_G * U(U'ΓU)^(-1)d
  # ============================================================================

  delta <- U %*% UtGammaU_inv %*% d
  b_PG <- b_RG + theta_G * as.vector(delta)

  # ============================================================================
  # COMPUTE METRICS
  # ============================================================================

  metrics <- .genomic_index_metrics(b_PG, Gamma, w, k_I, L_G,
    GAY = if (missing(GAY)) NULL else GAY,
    gmat = gmat
  )

  # Check proportionality: E should be proportional to d
  # For constrained traits (columns of U), compute U'E
  constrained_gains <- as.vector(t(U) %*% metrics$E_vec)
  gain_ratios <- constrained_gains / d
  gain_ratios[abs(d) < 1e-10] <- NA_real_ # Avoid division by zero

  # ============================================================================
  # BUILD SUMMARY OUTPUT
  # ============================================================================

  b_vec <- as.numeric(b_PG)
  b_vec <- round(b_vec, 4)

  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_along(b_vec)) # seq_len(length(b_vec)))

  summary_df <- data.frame(
    b_df,
    R = round(metrics$R, 4),
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    rHI = round(metrics$rHI, 4),
    theta_G = round(theta_G, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Name vectors
  trait_names <- colnames(Gamma)
  if (!is.null(trait_names)) {
    names(metrics$E_vec) <- trait_names
  }

  rownames(summary_df) <- NULL

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  result <- list(
    summary = summary_df,
    b = b_vec,
    E = metrics$E_vec,
    R = metrics$R,
    GA = metrics$GA,
    PRE = metrics$PRE,
    rHI = metrics$rHI,
    theta_G = theta_G,
    gain_ratios = gain_ratios,
    constrained_gains = constrained_gains,
    desired_gains = d,
    U = U,
    k_I = k_I,
    L_G = L_G
  )

  class(result) <- c("ppg_lgsi", "genomic_index", "list")

  result
}

# ==============================================================================
# 6.3 COMBINED RESTRICTED LINEAR GENOMIC SELECTION INDEX (CRLGSI)
# ==============================================================================

#' Combined Restricted Linear Genomic Selection Index (CRLGSI)
#'
#' @description
#' Implements the CRLGSI which combines phenotypic and genomic information
#' with restrictions on genetic gains. This extends CLGSI to include constraints.
#'
#' @param T_C Combined variance-covariance matrix (2t x 2t) where t = n_traits.
#'   Structure: [P, P_yg; P_yg', P_g] where P = phenotypic var, P_g = GEBV var,
#'   P_yg = covariance between phenotypes and GEBVs.
#'   Can be computed automatically if phen_mat and gebv_mat are provided.
#' @param Psi_C Combined genetic covariance matrix (2t x t).
#'   Structure: [G; C_gebv_g] where G = genetic var, C_gebv_g = Cov(GEBV, g).
#'   Can be computed automatically if gmat and reliability are provided.
#' @param phen_mat Optional. Matrix of phenotypes (n_genotypes x n_traits)
#' @param gebv_mat Optional. Matrix of GEBVs (n_genotypes x n_traits)
#' @param pmat Optional. Phenotypic variance-covariance matrix
#' @param gmat Optional. Genotypic variance-covariance matrix
#' @param wmat Economic weights matrix (n_traits x k), or vector
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param restricted_traits Vector of trait indices to restrict (default: NULL)
#' @param U Constraint matrix (2t x n_constraints for combined traits).
#'   Alternative to restricted_traits. Ignored if restricted_traits is provided.
#' @param reliability Optional. Reliability of GEBVs (r^2)
#' @param k_I Selection intensity (default: 2.063)
#' @param L_I Standardization constant (default: 1)
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients and metrics
#'     \item \code{b} - Vector of CRLGSI coefficients (β_CR)
#'     \item \code{b_y} - Coefficients for phenotypes
#'     \item \code{b_g} - Coefficients for GEBVs
#'     \item \code{E} - Expected genetic gains per trait
#'     \item \code{R} - Overall selection response
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 6, Section 6.3):}
#'
#' The CRLGSI combines phenotypic and genomic data with restrictions.
#'
#' Coefficient vector: \eqn{beta_CR = K_C * beta_C}
#'
#' Where K_C incorporates the restriction matrix.
#'
#' Selection response: \eqn{R_CR = (k_I / L_I) * sqrt(beta_CR' * T_C * beta_CR)}
#'
#' Expected gains: \eqn{E_CR = (k_I / L_I) * (Psi_C * beta_CR) / sqrt(beta_CR' * T_C * beta_CR)}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- 5
#'
#' phen_mat <- matrix(rnorm(n_genotypes * n_traits, 15, 3), n_genotypes, n_traits)
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits, 10, 2), n_genotypes, n_traits)
#'
#' gmat <- cov(phen_mat) * 0.6 # Genotypic component
#' pmat <- cov(phen_mat)
#'
#' w <- c(10, 8, 6, 4, 2)
#'
#' # Restrict traits 2 and 4
#' result <- crlgsi(
#'   phen_mat = phen_mat, gebv_mat = gebv_mat,
#'   pmat = pmat, gmat = gmat, wmat = w,
#'   restricted_traits = c(2, 4), reliability = 0.7
#' )
#' print(result$summary)
#' }
crlgsi <- function(T_C = NULL, Psi_C = NULL,
                   phen_mat = NULL, gebv_mat = NULL,
                   pmat = NULL, gmat = NULL,
                   wmat, wcol = 1,
                   restricted_traits = NULL, U = NULL,
                   reliability = NULL,
                   k_I = 2.063, L_I = 1, GAY = NULL) {
  # ============================================================================
  # INPUT VALIDATION AND MATRIX CONSTRUCTION
  # ============================================================================

  # Check if matrices provided directly or need to be computed
  has_direct_matrices <- !is.null(T_C) && !is.null(Psi_C)
  has_raw_data <- !is.null(phen_mat) && !is.null(gebv_mat)

  if (!has_direct_matrices && !has_raw_data) {
    stop("Must provide either (T_C, Psi_C) or (phen_mat, gebv_mat, pmat, gmat)")
  }

  # If raw data provided, compute T_C and Psi_C
  if (has_raw_data) {
    phen_mat <- as.matrix(phen_mat)
    gebv_mat <- as.matrix(gebv_mat)
    n_traits <- ncol(phen_mat)

    if (is.null(pmat)) pmat <- cov(phen_mat)
    if (is.null(gmat)) stop("gmat required when using raw data")

    pmat <- as.matrix(pmat)
    gmat <- as.matrix(gmat)

    # Compute covariance matrices
    P_gebv <- cov(gebv_mat)
    P_yg <- cov(phen_mat, gebv_mat)

    # Build T_C: [P, P_yg; P_yg', P_g]
    T_C <- rbind(
      cbind(pmat, P_yg),
      cbind(t(P_yg), P_gebv)
    )

    # Handle reliability for Psi_C
    if (is.null(reliability)) {
      reliability <- pmin(diag(P_gebv) / diag(gmat), 1)
    } else if (length(reliability) == 1) {
      reliability <- rep(reliability, n_traits)
    }

    # Build Psi_C: [G, Γ; Γ, Γ] (2t x 2t matrix)
    # Where G = genetic covariance, Γ = Cov(GEBV, g) = G * sqrt(reliability)
    Gamma_gebv_g <- sweep(gmat, 1, sqrt(reliability), "*")
    Psi_C <- rbind(
      cbind(gmat, Gamma_gebv_g),
      cbind(Gamma_gebv_g, Gamma_gebv_g)
    )
  } else {
    T_C <- as.matrix(T_C)
    Psi_C <- as.matrix(Psi_C)
    n_traits <- ncol(Psi_C) / 2 # Psi_C is 2t x 2t
  }

  if (nrow(T_C) != 2 * n_traits || ncol(T_C) != 2 * n_traits) {
    stop("T_C must be (2*n_traits x 2*n_traits)")
  }

  if (nrow(Psi_C) != 2 * n_traits || ncol(Psi_C) != 2 * n_traits) {
    stop("Psi_C must be (2*n_traits x 2*n_traits)")
  }

  # Handle weights
  wmat <- as.matrix(wmat)
  if (nrow(wmat) != n_traits) {
    stop("wmat must have ", n_traits, " rows")
  }
  w <- cpp_extract_vector(wmat, seq_len(n_traits), wcol - 1L)

  # Handle constraints
  if (!is.null(restricted_traits)) {
    if (!is.numeric(restricted_traits) || any(restricted_traits < 1) || any(restricted_traits > n_traits)) {
      stop("restricted_traits must be valid trait indices (1 to ", n_traits, ")")
    }
    # For combined index, must impose 2r restrictions (r traits + r GEBVs)
    # Build 2t x 2r constraint matrix
    r <- length(restricted_traits)
    U <- matrix(0, nrow = 2 * n_traits, ncol = 2 * r)
    for (i in seq_along(restricted_traits)) {
      trait_idx <- restricted_traits[i]
      # Constrain phenotype contribution (row trait_idx)
      U[trait_idx, i] <- 1
      # Constrain GEBV contribution (row n_traits + trait_idx)
      U[n_traits + trait_idx, r + i] <- 1
    }
  } else if (is.null(U)) {
    stop("Either 'restricted_traits' or 'U' must be provided")
  } else {
    U <- as.matrix(U)
  }

  # Constraint matrix for augmented system
  # The constraints are: U'*b_CR = 0 in the Psi_C metric
  # We solve [T_C, T_C*U; U'*T_C, 0][b; v] = [Psi_C*[w;w]; 0]

  U_TC <- U # Constraint matrix in combined variable space

  # ============================================================================
  # SOLVE AUGMENTED SYSTEM
  # ============================================================================

  n_constraints <- ncol(U_TC)

  # Build augmented matrix: [T_C, T_C*U_TC; U_TC'*T_C, 0]
  TC_UTC <- T_C %*% U_TC
  UTC_TC <- t(U_TC) %*% T_C
  zeros <- matrix(0, n_constraints, n_constraints)

  augmented_mat <- rbind(
    cbind(T_C, TC_UTC),
    cbind(UTC_TC, zeros)
  )

  # Build RHS: [Psi_C*a_C; 0] where a_C = [w; 0]
  # Economic weights apply only to breeding values, not GEBVs
  a_C <- c(w, rep(0, n_traits))
  Psi_w <- Psi_C %*% a_C
  rhs <- c(Psi_w, rep(0, n_constraints))

  # Solve using ginv for robustness
  solution <- ginv(augmented_mat) %*% rhs

  # Extract β_CR
  b_CR <- solution[1:(2 * n_traits)]

  if (any(is.na(b_CR)) || any(is.infinite(b_CR))) {
    stop("CRLGSI coefficients contain NA or Inf. Check matrix conditioning.")
  }

  # Split into phenotype and GEBV coefficients
  b_y <- b_CR[1:n_traits]
  b_g <- b_CR[(n_traits + 1):(2 * n_traits)]

  # ============================================================================
  # COMPUTE METRICS
  # ============================================================================

  metrics <- .combined_index_metrics(b_CR, T_C, Psi_C, w, k_I, L_I,
    GAY = if (missing(GAY)) NULL else GAY
  )

  # Check constraint satisfaction (constraints are on both phenotype and GEBV)
  # Extract genetic gains for the r restricted traits
  if (!is.null(restricted_traits)) {
    constrained_response <- metrics$E_vec[restricted_traits]
  } else {
    # For custom U, need to project back to trait space
    # This is more complex, for now just report NA
    constrained_response <- rep(NA_real_, ncol(U) / 2)
  }

  # ============================================================================
  # BUILD SUMMARY OUTPUT
  # ============================================================================

  b_y_vec <- round(as.numeric(b_y), 4)
  b_g_vec <- round(as.numeric(b_g), 4)

  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_along(b_y_vec)) # seq_len(length(b_y_vec)))

  b_g_df <- as.data.frame(matrix(b_g_vec, nrow = 1))
  colnames(b_g_df) <- paste0("b_g.", seq_along(b_g_vec)) # seq_len(length(b_g_vec)))

  summary_df <- data.frame(
    b_y_df,
    b_g_df,
    R = round(metrics$R, 4),
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    rHI = round(metrics$rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(summary_df) <- NULL

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  result <- list(
    summary = summary_df,
    b = as.numeric(b_CR),
    b_y = b_y_vec,
    b_g = b_g_vec,
    E = metrics$E_vec,
    R = metrics$R,
    GA = metrics$GA,
    PRE = metrics$PRE,
    rHI = metrics$rHI,
    sigma_I = metrics$sigma_I,
    U = U,
    constrained_response = constrained_response,
    k_I = k_I,
    L_I = L_I
  )

  class(result) <- c("crlgsi", "genomic_index", "list")

  result
}

# ==============================================================================
# 6.4 COMBINED PREDETERMINED PROPORTIONAL GAINS LGSI (CPPG-LGSI)
# ==============================================================================

#' Combined Predetermined Proportional Gains Linear Genomic Selection Index (CPPG-LGSI)
#'
#' @description
#' Implements the CPPG-LGSI which combines phenotypic and genomic information
#' while achieving predetermined proportional gains between traits.
#' This is the most general constrained genomic index.
#'
#' @param T_C Combined variance-covariance matrix (2t x 2t)
#' @param Psi_C Combined genetic covariance matrix (2t x t)
#' @param d Vector of desired proportional gains (length n_traits or n_constraints)
#' @param phen_mat Optional. Matrix of phenotypes for automatic T_C computation
#' @param gebv_mat Optional. Matrix of GEBVs for automatic T_C computation
#' @param pmat Optional. Phenotypic variance-covariance matrix
#' @param gmat Optional. Genotypic variance-covariance matrix
#' @param wmat Optional. Economic weights for GA/PRE calculation
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param U Optional. Constraint matrix (n_traits x n_constraints)
#' @param reliability Optional. Reliability of GEBVs (r^2)
#' @param k_I Selection intensity (default: 2.063)
#' @param L_I Standardization constant (default: 1)
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summary} - Data frame with coefficients and metrics
#'     \item \code{b} - Vector of CPPG-LGSI coefficients (β_CP)
#'     \item \code{b_y} - Coefficients for phenotypes
#'     \item \code{b_g} - Coefficients for GEBVs
#'     \item \code{E} - Expected genetic gains per trait
#'     \item \code{theta_CP} - Proportionality constant
#'     \item \code{gain_ratios} - Ratios of achieved to desired gains
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 6, Section 6.4):}
#'
#' Coefficient vector: \eqn{beta_CP = beta_CR + theta_CP * delta_CP}
#'
#' Where beta_CR is from CRLGSI and:
#'
#' \deqn{theta_CP = (beta_C' * Phi_C * (Phi_C' * T_C^{-1} * Phi_C)^{-1} * d_C) / (d_C' * (Phi_C' * T_C^{-1} * Phi_C)^{-1} * d_C)}
#'
#' Selection response: \eqn{R_CP = (k_I / L_I) * sqrt(beta_CP' * T_C * beta_CP)}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- 5
#'
#' phen_mat <- matrix(rnorm(n_genotypes * n_traits, 15, 3), n_genotypes, n_traits)
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits, 10, 2), n_genotypes, n_traits)
#'
#' gmat <- cov(phen_mat) * 0.6
#' pmat <- cov(phen_mat)
#'
#' # Desired proportional gains
#' d <- c(2, 1, 1, 0.5, 0)
#'
#' w <- c(10, 8, 6, 4, 2)
#'
#' result <- cppg_lgsi(
#'   phen_mat = phen_mat, gebv_mat = gebv_mat,
#'   pmat = pmat, gmat = gmat, d = d, wmat = w,
#'   reliability = 0.7
#' )
#' print(result$summary)
#' print(result$gain_ratios)
#' }
cppg_lgsi <- function(T_C = NULL, Psi_C = NULL, d,
                      phen_mat = NULL, gebv_mat = NULL,
                      pmat = NULL, gmat = NULL,
                      wmat = NULL, wcol = 1, U = NULL,
                      reliability = NULL,
                      k_I = 2.063, L_I = 1, GAY = NULL) {
  # ============================================================================
  # INPUT VALIDATION AND MATRIX CONSTRUCTION
  # ============================================================================

  d <- as.numeric(d)

  # Check if matrices provided directly or need to be computed
  has_direct_matrices <- !is.null(T_C) && !is.null(Psi_C)
  has_raw_data <- !is.null(phen_mat) && !is.null(gebv_mat)

  if (!has_direct_matrices && !has_raw_data) {
    stop("Must provide either (T_C, Psi_C) or (phen_mat, gebv_mat, pmat, gmat)")
  }

  # Compute T_C and Psi_C if needed
  if (has_raw_data) {
    phen_mat <- as.matrix(phen_mat)
    gebv_mat <- as.matrix(gebv_mat)
    n_traits <- ncol(phen_mat)

    if (is.null(pmat)) pmat <- cov(phen_mat)
    if (is.null(gmat)) stop("gmat required when using raw data")

    pmat <- as.matrix(pmat)
    gmat <- as.matrix(gmat)

    P_gebv <- cov(gebv_mat)
    P_yg <- cov(phen_mat, gebv_mat)

    T_C <- rbind(
      cbind(pmat, P_yg),
      cbind(t(P_yg), P_gebv)
    )

    if (is.null(reliability)) {
      reliability <- pmin(diag(P_gebv) / diag(gmat), 1)
    } else if (length(reliability) == 1) {
      reliability <- rep(reliability, n_traits)
    }

    C_gebv_g <- sweep(gmat, 1, sqrt(reliability), "*")
    Psi_C <- rbind(
      cbind(gmat, C_gebv_g),
      cbind(C_gebv_g, C_gebv_g)
    )
  } else {
    T_C <- as.matrix(T_C)
    Psi_C <- as.matrix(Psi_C)
    n_traits <- ncol(Psi_C) / 2 # Psi_C is 2t x 2t
  }

  # Handle U matrix and double d for combined indices
  if (is.null(U)) {
    if (length(d) != n_traits) {
      stop("d must have length ", n_traits, " when U is not provided")
    }
    # For combined indices, d_C = [d; d] and U must impose 2r constraints
    # Build 2t x 2t identity-like constraint matrix
    U_base <- diag(n_traits)
    # Stack: [U_base in top t rows, U_base in bottom t rows]
    # This gives 2t x t matrix, but we need 2t x 2t for 2r constraints
    U <- rbind(
      cbind(U_base, matrix(0, n_traits, n_traits)),
      cbind(matrix(0, n_traits, n_traits), U_base)
    )
    # Double the d vector: d_C = [d; d]
    d_C <- c(d, d)
  } else {
    U <- as.matrix(U)
    if (nrow(U) != 2 * n_traits) {
      stop("U must have ", 2 * n_traits, " rows for combined indices")
    }
    # U should have 2r columns for r constraints
    # d can be provided as r-length (gets doubled) or 2r-length (used as-is)
    r <- ncol(U) / 2
    if (length(d) == r) {
      # Auto-double if user provided single d for r traits
      d_C <- c(d, d)
    } else if (length(d) == ncol(U)) {
      # User provided full d_C
      d_C <- d
    } else {
      stop("d must have length ", r, " or ", ncol(U), " when providing custom U")
    }
  }

  # Handle weights
  w <- NULL
  if (!is.null(wmat)) {
    wmat <- as.matrix(wmat)
    if (nrow(wmat) != n_traits) {
      stop("wmat must have ", n_traits, " rows")
    }
    w <- cpp_extract_vector(wmat, seq_len(n_traits), wcol - 1L)
  } else {
    w <- rep(1, n_traits) # Default equal weights
  }

  # ============================================================================
  # STEP 1: Compute β_CR (CRLGSI component)
  # ============================================================================

  # For CPPG-LGSI, first solve CRLGSI with all traits restricted
  # Build constraint matrix for all traits (2r = 2t constraints)
  all_traits <- seq_len(n_traits)
  U_crlgsi <- matrix(0, nrow = 2 * n_traits, ncol = 2 * n_traits)
  for (i in seq_along(all_traits)) {
    trait_idx <- all_traits[i]
    U_crlgsi[trait_idx, i] <- 1
    U_crlgsi[n_traits + trait_idx, n_traits + i] <- 1
  }

  crlgsi_result <- crlgsi(
    T_C = T_C, Psi_C = Psi_C, wmat = w, wcol = 1,
    U = U_crlgsi, k_I = k_I, L_I = L_I, GAY = NULL
  )
  b_CR <- crlgsi_result$b

  # ============================================================================
  # STEP 2: Compute Proportionality Constant θ_CP
  # ============================================================================

  # Φ_C = Psi_C * U (maps constraint space to combined variable space)
  Phi_C <- Psi_C %*% U

  # T_C^(-1) using ginv
  T_C_inv <- ginv(T_C)

  # (Φ_C' T_C^(-1) Φ_C)
  Phi_TC_Phi <- t(Phi_C) %*% T_C_inv %*% Phi_C

  # (Φ_C' T_C^(-1) Φ_C)^(-1) using ginv for robustness
  Phi_TC_Phi_inv <- ginv(Phi_TC_Phi)

  # β_C = T_C^(-1) Psi_C a_C where a_C = [w; 0] (unrestricted combined index)
  # Economic weights apply only to breeding values, not GEBVs
  a_C <- c(w, rep(0, n_traits))
  beta_C <- T_C_inv %*% Psi_C %*% a_C

  # θ_CP numerator: β_C' Φ_C (Φ_C' T_C^(-1) Φ_C)^(-1) d_C
  theta_num <- as.numeric(t(beta_C) %*% Phi_C %*% Phi_TC_Phi_inv %*% d_C)

  # θ_CP denominator: d_C' (Φ_C' T_C^(-1) Φ_C)^(-1) d_C
  theta_denom <- as.numeric(t(d_C) %*% Phi_TC_Phi_inv %*% d_C)

  theta_CP <- if (abs(theta_denom) > 1e-10) theta_num / theta_denom else 0

  # ============================================================================
  # STEP 3: Compute β_CP = β_CR + θ_CP * δ_CP
  # ============================================================================

  # δ_CP = T_C^(-1) Φ_C (Φ_C' T_C^(-1) Φ_C)^(-1) d_C
  delta_CP <- T_C_inv %*% Phi_C %*% Phi_TC_Phi_inv %*% d_C

  b_CP <- b_CR + theta_CP * as.vector(delta_CP)

  # Split into phenotype and GEBV coefficients
  b_y <- b_CP[1:n_traits]
  b_g <- b_CP[(n_traits + 1):(2 * n_traits)]

  # ============================================================================
  # STEP 4: Compute Metrics
  # ============================================================================

  metrics <- .combined_index_metrics(b_CP, T_C, Psi_C, w, k_I, L_I,
    GAY = if (missing(GAY)) NULL else GAY
  )

  # Check proportionality (extract first t elements of gains)
  constrained_gains <- metrics$E_vec[1:n_traits]
  if (length(d) == n_traits) {
    gain_ratios <- constrained_gains / d
    gain_ratios[abs(d) < 1e-10] <- NA_real_
  } else {
    # Custom U provided with d of length r (not n_traits) — ratio not comparable
    gain_ratios <- rep(NA_real_, n_traits)
  }

  # ============================================================================
  # BUILD SUMMARY OUTPUT
  # ============================================================================

  b_y_vec <- round(as.numeric(b_y), 4)
  b_g_vec <- round(as.numeric(b_g), 4)

  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_along(b_y_vec)) # seq_len(length(b_y_vec)))

  b_g_df <- as.data.frame(matrix(b_g_vec, nrow = 1))
  colnames(b_g_df) <- paste0("b_g.", seq_along(b_g_vec)) # seq_len(length(b_g_vec)))

  summary_df <- data.frame(
    b_y_df,
    b_g_df,
    R = round(metrics$R, 4),
    GA = round(metrics$GA, 4),
    PRE = round(metrics$PRE, 4),
    rHI = round(metrics$rHI, 4),
    theta_CP = round(theta_CP, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(summary_df) <- NULL

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  result <- list(
    summary = summary_df,
    b = as.numeric(b_CP),
    b_y = b_y_vec,
    b_g = b_g_vec,
    E = metrics$E_vec,
    R = metrics$R,
    GA = metrics$GA,
    PRE = metrics$PRE,
    rHI = metrics$rHI,
    theta_CP = theta_CP,
    gain_ratios = gain_ratios,
    constrained_gains = constrained_gains,
    desired_gains = d,
    U = U,
    k_I = k_I,
    L_I = L_I
  )

  class(result) <- c("cppg_lgsi", "genomic_index", "list")

  result
}
