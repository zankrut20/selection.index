#' Multistage Linear Genomic Selection Indices (Chapter 9)
#' @name multistage_genomic_indices
#'
#' @description
#' Implements multistage genomic selection index methods from Chapter 9.
#' These methods combine genomic estimated breeding values (GEBVs) across
#' multiple stages with covariance adjustments due to selection effects
#' using Cochran/Cunningham's method.
#'
#' Methods included:
#' - MLGSI: Multistage Linear Genomic Selection Index (Section 9.4)
#' - MRLGSI: Multistage Restricted Linear Genomic Selection Index (Section 9.5)
#' - MPPG-LGSI: Multistage Predetermined Proportional Gain LGSI (Section 9.6)
#'
#' All implementations use C++ primitives for mathematical operations.
#'
#' @references
#' Cochran, W. G. (1951). Improvement by means of selection.
#' Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability, 449-470.
#'
#' Young, S. S. Y. (1964). Multi-stage selection for genetic gain.
#' Heredity, 19(1), 131-144.
#'
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 9.
#'
#' @keywords internal
#' @importFrom stats setNames pnorm qnorm
#' @importFrom MASS ginv
NULL

# ==============================================================================
# HELPER FUNCTIONS (Reuse from multistage_phenotypic_indices.R)
# ==============================================================================

#' Compute Cochran/Cunningham covariance adjustment for genomic matrices
#' @keywords internal
#' @noRd
.cochran_adjustment_genomic <- function(Gamma, Gamma1, beta1, A, k1, tau) {
  # Gamma* = Gamma - u * A beta1 beta1' A' / (beta1'Gamma1 beta1)
  # where A is n x n1 (covariance between all n GEBVs and n1 stage 1 true BVs)
  u <- k1 * (k1 - tau)
  beta1_Gamma1_beta1 <- cpp_quadratic_form_sym(beta1, Gamma1)

  if (beta1_Gamma1_beta1 <= 0) {
    warning("Invalid genomic variance at stage 1. Returning unadjusted matrix.")
    return(Gamma)
  }

  # A is n x n1, beta1 is n1 x 1
  # A %*% beta1 is n x 1
  # crossprod(beta1, t(A)) = t(beta1) %*% t(A) = t(A %*% beta1) is 1 x n
  # (n x 1) %*% (1 x n) = n x n
  A_beta1 <- A %*% beta1 # n x 1
  beta1_A <- crossprod(beta1, t(A)) # 1 x n
  adjustment <- u * (A_beta1 %*% beta1_A) / beta1_Gamma1_beta1
  Gamma_star <- Gamma - adjustment

  Gamma_star
}

#' Compute genomic stage metrics for multistage indices
#' @keywords internal
#' @noRd
.genomic_stage_metrics <- function(beta, Gamma, A, w, k) {
  beta <- as.numeric(beta)
  beta_Gamma_beta <- cpp_quadratic_form_sym(beta, Gamma)

  # Handle NaN/NA values
  if (!is.finite(beta_Gamma_beta)) {
    warning("Invalid genomic variance (NaN or NA). Returning NA metrics.")
    return(list(
      sigma_I = NA_real_,
      R = NA_real_,
      E = rep(NA_real_, if (!is.null(A)) nrow(A) else length(w)),
      beta_Gamma_beta = beta_Gamma_beta
    ))
  }

  sigma_I <- if (beta_Gamma_beta > 0) sqrt(beta_Gamma_beta) else NA_real_

  # Selection response: R = k * sqrt(beta'Gamma beta)
  R <- if (!is.na(sigma_I)) k * sigma_I else NA_real_

  # Expected genetic gain per trait: E = k * A'beta / sqrt(beta'Gamma beta)
  # Or for stage 2 with w: E = k * Gamma*w / sqrt(w'Gamma*w)
  if (!is.null(A)) {
    A_beta <- A %*% beta
    E <- if (!is.na(sigma_I) && sigma_I > 0) k * A_beta / sigma_I else rep(NA_real_, length(A_beta))
  } else {
    # Stage 2 case where beta = w
    Gamma_w <- Gamma %*% w
    E <- if (!is.na(sigma_I) && sigma_I > 0) k * Gamma_w / sigma_I else rep(NA_real_, length(Gamma_w))
  }

  # Accuracy: rho_HI = sqrt(beta'Gamma beta / w'C w)
  # This will be computed separately with C matrix

  list(
    sigma_I = sigma_I,
    R = R,
    E = as.numeric(E),
    beta_Gamma_beta = beta_Gamma_beta
  )
}

#' Compute genomic index correlation
#' @keywords internal
#' @noRd
.genomic_index_correlation <- function(beta1, beta2, Gamma1, Gamma, A) {
  # rho_I1I2 = beta1' A' beta2 / sqrt(beta1'Gamma1 beta1) sqrt(beta2'Gamma beta2)
  # For MLGSI, beta2 = w (economic weights for all traits)
  # A is n x n1 (covariance between all n stage 2 GEBVs and n1 stage 1 true BVs)
  # So A' is n1 x n, beta1 is n1 x 1, beta2 is n x 1
  # beta1' A' beta2 = (1 x n1) * (n1 x n) * (n x 1) = scalar

  beta1_Gamma1_beta1 <- cpp_quadratic_form_sym(beta1, Gamma1)
  beta2_Gamma_beta2 <- cpp_quadratic_form_sym(beta2, Gamma)

  if (beta1_Gamma1_beta1 <= 0 || beta2_Gamma_beta2 <= 0) {
    warning("Invalid variance for genomic correlation calculation.")
    return(NA_real_)
  }


  numerator <- as.numeric(crossprod(beta1, crossprod(A, beta2)))
  denominator <- sqrt(beta1_Gamma1_beta1) * sqrt(beta2_Gamma_beta2)

  rho <- numerator / denominator
  rho
}

#' Young's method for selection intensities (same as phenotypic)
#' @keywords internal
#' @noRd
.young_intensities <- function(p, rho_12) {
  if (p <= 0 || p >= 1) {
    stop("Selection proportion p must be between 0 and 1")
  }

  # Handle edge cases for correlation
  rho_12 <- max(min(rho_12, 0.999), -0.999) # Clamp to valid range

  c1 <- qnorm(1 - p)
  c3 <- qnorm(1 - sqrt(p))

  z <- function(x) dnorm(x)
  Q <- function(x) 1 - pnorm(x)

  a <- (1 - rho_12) / sqrt(2 * (1 - rho_12))
  b <- (1 + rho_12) / sqrt(2 * (1 + rho_12))

  k1 <- (z(c1) * Q(a) / p) + (z(c3) * Q(b) * sqrt((1 + rho_12) / 2) / p)
  k2 <- (rho_12 * z(c1) * Q(a) / p) + (z(c3) * Q(b) * sqrt((1 + rho_12) / 2) / p)

  list(k1 = k1, k2 = k2)
}

# ==============================================================================
# 9.4 MULTISTAGE LINEAR GENOMIC SELECTION INDEX (MLGSI)
# ==============================================================================

#' Multistage Linear Genomic Selection Index (MLGSI)
#'
#' @description
#' Implements the two-stage Linear Genomic Selection Index where selection
#' is based on GEBVs at both stages with covariance adjustments due to
#' selection effects.
#'
#' @param Gamma1 GEBV variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param Gamma GEBV variance-covariance matrix for all traits at stage 2 (n x n)
#' @param A1 Covariance matrix between GEBVs and true breeding values for stage 1 (n1 x n1)
#' @param A Covariance matrix between GEBVs and true breeding values for stage 2 (n x n1)
#' @param C Genotypic variance-covariance matrix for all traits (n x n)
#' @param G1 Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param P1 Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param wmat Economic weights vector or matrix (n x k)
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param selection_proportion Proportion selected at each stage (default: 0.1)
#' @param use_young_method Logical. Use Young's method for selection intensities (default: FALSE).
#'   Young's method tends to overestimate intensities; manual intensities are recommended.
#' @param k1_manual Manual selection intensity for stage 1
#' @param k2_manual Manual selection intensity for stage 2
#' @param tau Standardized truncation point
#' @param reliability Optional reliability vector for computing A matrices
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{beta1} - Stage 1 genomic index coefficients
#'     \item \code{w} - Economic weights (stage 2 coefficients)
#'     \item \code{stage1_metrics} - List with stage 1 metrics (R1, E1, rho_HI1)
#'     \item \code{stage2_metrics} - List with stage 2 metrics (R2, E2, rho_HI2)
#'     \item \code{Gamma_star} - Adjusted genomic covariance matrix at stage 2
#'     \item \code{C_star} - Adjusted genotypic covariance matrix at stage 2
#'     \item \code{rho_I1I2} - Correlation between stage 1 and stage 2 indices
#'     \item \code{k1} - Selection intensity at stage 1
#'     \item \code{k2} - Selection intensity at stage 2
#'     \item \code{summary_stage1} - Data frame with stage 1 summary
#'     \item \code{summary_stage2} - Data frame with stage 2 summary
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' Stage 1: The genomic index is \eqn{I_1 = \mathbf{\beta}_1' \mathbf{\gamma}_1}
#'
#' Coefficients: \eqn{\mathbf{\beta}_1 = \mathbf{\Gamma}_1^{-1}\mathbf{A}_1\mathbf{w}_1}
#'
#' Stage 2: The index uses economic weights directly: \eqn{I_2 = \mathbf{w}' \mathbf{\gamma}}
#'
#' Adjusted genomic covariance matrix:
#' \deqn{\mathbf{\Gamma}^* = \mathbf{\Gamma} - u \frac{\mathbf{A}_1'\mathbf{\beta}_1\mathbf{\beta}_1'\mathbf{A}_1}{\mathbf{\beta}_1'\mathbf{\Gamma}_1\mathbf{\beta}_1}}
#'
#' Adjusted genotypic covariance matrix:
#' \deqn{\mathbf{C}^* = \mathbf{C} - u \frac{\mathbf{G}_1'\mathbf{b}_1\mathbf{b}_1'\mathbf{G}_1}{\mathbf{b}_1'\mathbf{P}_1\mathbf{b}_1}}
#'
#' where \eqn{u = k_1(k_1 - \tau)}
#'
#' Accuracy at stage 1: \eqn{\rho_{HI_1} = \sqrt{\frac{\mathbf{\beta}_1'\mathbf{\Gamma}_1\mathbf{\beta}_1}{\mathbf{w}'\mathbf{C}\mathbf{w}}}}
#'
#' Accuracy at stage 2: \eqn{\rho_{HI_2} = \sqrt{\frac{\mathbf{w}'\mathbf{\Gamma}^*\mathbf{w}}{\mathbf{w}'\mathbf{C}^*\mathbf{w}}}}
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 9, Section 9.4.
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-stage genomic selection example
#' # Stage 1: Select based on GEBVs for 3 traits
#' # Stage 2: Select based on GEBVs for all 7 traits
#'
#' # Compute covariance matrices
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate GEBV covariances (in practice, compute from genomic prediction)
#' set.seed(123)
#' reliability <- 0.7
#' Gamma1 <- reliability * gmat[1:3, 1:3]
#' Gamma <- reliability * gmat
#' A1 <- reliability * gmat[1:3, 1:3]
#' A <- gmat[, 1:3]
#'
#' # Economic weights
#' weights <- c(10, 8, 6, 4, 3, 2, 1)
#'
#' # Run MLGSI
#' result <- mlgsi(
#'   Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
#'   C = gmat, G1 = gmat[1:3, 1:3], P1 = pmat[1:3, 1:3],
#'   wmat = weights, selection_proportion = 0.1
#' )
#'
#' print(result$summary_stage1)
#' print(result$summary_stage2)
#' }
mlgsi <- function(Gamma1, Gamma, A1, A, C, G1, P1, wmat, wcol = 1,
                  selection_proportion = 0.1,
                  use_young_method = FALSE,
                  k1_manual = 2.063,
                  k2_manual = 2.063,
                  tau = NULL,
                  reliability = NULL) {
  # ============================================================================
  # STEP 1: Input Validation
  # ============================================================================

  Gamma1 <- as.matrix(Gamma1)
  Gamma <- as.matrix(Gamma)
  A1 <- as.matrix(A1)
  A <- as.matrix(A)
  C <- as.matrix(C)
  G1 <- as.matrix(G1)
  P1 <- as.matrix(P1)

  n1 <- nrow(Gamma1)
  n <- nrow(Gamma)

  # Process weights
  wmat <- as.matrix(wmat)
  if (ncol(wmat) == 1) {
    w <- as.numeric(wmat)
  } else {
    w <- as.numeric(wmat[, wcol])
  }

  if (length(w) != n) {
    stop("Weight vector length must match number of traits at stage 2")
  }

  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  # ============================================================================
  # STEP 2: Stage 1 - Compute genomic index coefficients
  # ============================================================================

  # beta1 = Gamma1^{-1} A1 w1
  w1 <- w[1:n1]
  Gamma1_inv_A1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    Gamma1_inv_A1[, j] <- cpp_symmetric_solve(Gamma1, A1[, j])
  }
  beta1 <- Gamma1_inv_A1 %*% w1

  # ============================================================================
  # STEP 3: Stage 2 - Coefficients are just economic weights
  # ============================================================================

  # For MLGSI, stage 2 uses w directly
  beta2 <- w

  # ============================================================================
  # STEP 4: Compute phenotypic index coefficients for C* adjustment
  # ============================================================================

  # For adjusting C*, we need phenotypic index from stage 1
  # b1 = P1^{-1} G1 w1
  P1_inv_G1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    P1_inv_G1[, j] <- cpp_symmetric_solve(P1, G1[, j])
  }
  b1 <- P1_inv_G1 %*% w1

  # ============================================================================
  # STEP 5: Compute correlation between genomic indices
  # ============================================================================

  rho_I1I2 <- .genomic_index_correlation(beta1, beta2, Gamma1, Gamma, A)

  # ============================================================================
  # STEP 6: Compute selection intensities
  # ============================================================================

  # Set tau from selection proportion if not provided
  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  if (use_young_method && !is.na(rho_I1I2)) {
    intensities <- tryCatch(
      {
        .young_intensities(selection_proportion, rho_I1I2)
      },
      error = function(e) {
        warning("Young's method failed. Using manual intensities.")
        list(k1 = k1_manual, k2 = k2_manual)
      }
    )
    k1 <- intensities$k1
    k2 <- intensities$k2
  } else {
    k1 <- k1_manual
    k2 <- k2_manual
  }

  # ============================================================================
  # STEP 7: Adjust genomic and genotypic covariance matrices
  # ============================================================================

  # Adjust Gamma for stage 2
  Gamma_star <- .cochran_adjustment_genomic(Gamma, Gamma1, beta1, A, k1, tau)

  # Adjust C for stage 2 (using phenotypic index)
  u <- k1 * (k1 - tau)
  b1Pb1 <- cpp_quadratic_form_sym(b1, P1)

  if (b1Pb1 > 0 && nrow(b1) == ncol(G1)) {
    # C* = C - u * C[,1:n1] b1 b1' C[1:n1,] / (b1'P1b1)
    n1 <- nrow(G1)
    C_col_1_n1 <- C[, 1:n1]
    C_row_1_n1 <- C[1:n1, ]

    C_col_b1 <- C_col_1_n1 %*% b1 # n x 1
    b1_C_row <- crossprod(b1, C_row_1_n1) # 1 x n
    adjustment_C <- u * (C_col_b1 %*% b1_C_row) / b1Pb1 # n x n
    C_star <- C - adjustment_C
  } else {
    warning("Invalid phenotypic variance at stage 1. Using unadjusted C.")
    C_star <- C
  }

  # ============================================================================
  # STEP 8: Compute stage metrics
  # ============================================================================

  # Stage 1 metrics
  stage1_metrics <- .genomic_stage_metrics(beta1, Gamma1, A1, w1, k1)

  # Compute rho_HI1
  wCw <- cpp_quadratic_form_sym(w, C)
  rho_HI1 <- if (stage1_metrics$beta_Gamma_beta > 0 && wCw > 0) {
    sqrt(stage1_metrics$beta_Gamma_beta / wCw)
  } else {
    NA_real_
  }
  stage1_metrics$rho_HI <- rho_HI1

  # Stage 2 metrics (using w as coefficients)
  stage2_metrics <- .genomic_stage_metrics(w, Gamma_star, NULL, w, k2)

  # Compute rho_HI2
  wCstarw <- cpp_quadratic_form_sym(w, C_star)
  rho_HI2 <- if (stage2_metrics$beta_Gamma_beta > 0 && wCstarw > 0) {
    sqrt(stage2_metrics$beta_Gamma_beta / wCstarw)
  } else {
    NA_real_
  }
  stage2_metrics$rho_HI <- rho_HI2

  # ============================================================================
  # STEP 9: Create summary data frames
  # ============================================================================

  trait_names_1 <- if (!is.null(colnames(Gamma1))) colnames(Gamma1) else paste0("Trait", 1:n1)
  trait_names_2 <- if (!is.null(colnames(Gamma))) colnames(Gamma) else paste0("Trait", 1:n)

  summary_stage1 <- data.frame(
    Stage = 1,
    Trait = trait_names_1,
    beta = round(as.numeric(beta1), 4),
    E = round(stage1_metrics$E, 4),
    stringsAsFactors = FALSE
  )

  summary_stage2 <- data.frame(
    Stage = 2,
    Trait = trait_names_2,
    w = round(as.numeric(w), 4),
    E = round(stage2_metrics$E, 4),
    stringsAsFactors = FALSE
  )

  # ============================================================================
  # STEP 10: Return results
  # ============================================================================

  result <- list(
    beta1 = as.numeric(beta1),
    w = as.numeric(w),
    stage1_metrics = list(
      R = stage1_metrics$R,
      E = setNames(stage1_metrics$E, trait_names_1),
      rho_HI = stage1_metrics$rho_HI,
      sigma_I = stage1_metrics$sigma_I
    ),
    stage2_metrics = list(
      R = stage2_metrics$R,
      E = setNames(stage2_metrics$E, trait_names_2),
      rho_HI = stage2_metrics$rho_HI,
      sigma_I = stage2_metrics$sigma_I
    ),
    Gamma_star = Gamma_star,
    C_star = C_star,
    rho_I1I2 = rho_I1I2,
    k1 = k1,
    k2 = k2,
    tau = tau,
    selection_proportion = selection_proportion,
    summary_stage1 = summary_stage1,
    summary_stage2 = summary_stage2
  )

  class(result) <- c("mlgsi", "multistage_genomic_index", "list")

  result
}

# ==============================================================================
# 9.5 MULTISTAGE RESTRICTED LINEAR GENOMIC SELECTION INDEX (MRLGSI)
# ==============================================================================

#' Multistage Restricted Linear Genomic Selection Index (MRLGSI)
#'
#' @description
#' Implements the two-stage Restricted Linear Genomic Selection Index where
#' certain traits are constrained to have zero genetic gain at each stage
#' using GEBVs.
#'
#' @param Gamma1 GEBV variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param Gamma GEBV variance-covariance matrix for all traits at stage 2 (n x n)
#' @param A1 Covariance matrix between GEBVs and true breeding values for stage 1 (n1 x n1)
#' @param A Covariance matrix between GEBVs and true breeding values for stage 2 (n x n1)
#' @param C Genotypic variance-covariance matrix for all traits (n x n)
#' @param G1 Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param P1 Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param wmat Economic weights vector or matrix (n x k)
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param C1 Constraint matrix for stage 1 (n1 x r1)
#' @param C2 Constraint matrix for stage 2 (n x r2)
#' @param selection_proportion Proportion selected at each stage (default: 0.1)
#' @param use_young_method Logical. Use Young's method for selection intensities (default: FALSE).
#'   Young's method tends to overestimate intensities; manual intensities are recommended.
#' @param k1_manual Manual selection intensity for stage 1
#' @param k2_manual Manual selection intensity for stage 2
#' @param tau Standardized truncation point
#'
#' @return List with components similar to mlgsi, plus:
#'   \itemize{
#'     \item \code{beta_R1} - Restricted stage 1 coefficients
#'     \item \code{beta_R2} - Restricted stage 2 coefficients
#'     \item \code{K_G1} - Restriction matrix for stage 1
#'     \item \code{K_G2} - Restriction matrix for stage 2
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The restricted genomic coefficients are:
#' \deqn{\mathbf{\beta}_{R_1} = \mathbf{K}_{G_1}\mathbf{\beta}_1}
#' \deqn{\mathbf{\beta}_{R_2} = \mathbf{K}_{G_2}\mathbf{w}}
#'
#' where restriction matrices are computed similarly to RLGSI
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 9, Section 9.5.
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-stage restricted genomic selection
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' reliability <- 0.7
#' Gamma1 <- reliability * gmat[1:3, 1:3]
#' Gamma <- reliability * gmat
#' A1 <- reliability * gmat[1:3, 1:3]
#' A <- gmat[, 1:3]
#'
#' # Constraint matrices
#' C1 <- matrix(0, nrow = 3, ncol = 1)
#' C1[1, 1] <- 1 # Restrict trait 1 at stage 1
#'
#' C2 <- matrix(0, nrow = 7, ncol = 2)
#' C2[1, 1] <- 1 # Restrict trait 1 at stage 2
#' C2[3, 2] <- 1 # Restrict trait 3 at stage 2
#'
#' weights <- c(10, 8, 6, 4, 3, 2, 1)
#'
#' result <- mrlgsi(
#'   Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
#'   C = gmat, G1 = gmat[1:3, 1:3], P1 = pmat[1:3, 1:3],
#'   wmat = weights, C1 = C1, C2 = C2
#' )
#' }
mrlgsi <- function(Gamma1, Gamma, A1, A, C, G1, P1, wmat, wcol = 1,
                   C1, C2,
                   selection_proportion = 0.1,
                   use_young_method = FALSE,
                   k1_manual = 2.063,
                   k2_manual = 2.063,
                   tau = NULL) {
  # ============================================================================
  # STEP 1: Input Validation
  # ============================================================================

  Gamma1 <- as.matrix(Gamma1)
  Gamma <- as.matrix(Gamma)
  A1 <- as.matrix(A1)
  C1 <- as.matrix(C1)
  C2 <- as.matrix(C2)

  n1 <- nrow(Gamma1)
  n <- nrow(Gamma)

  # Process weights
  wmat <- as.matrix(wmat)
  if (ncol(wmat) == 1) {
    w <- as.numeric(wmat)
  } else {
    w <- as.numeric(wmat[, wcol])
  }

  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  # ============================================================================
  # STEP 2: Compute unrestricted genomic coefficients
  # ============================================================================

  w1 <- w[1:n1]

  # Stage 1: beta1 = Gamma1^{-1} A1 w1
  Gamma1_inv_A1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    Gamma1_inv_A1[, j] <- cpp_symmetric_solve(Gamma1, A1[, j])
  }
  beta1 <- Gamma1_inv_A1 %*% w1

  # Stage 2: beta2 = w (for MLGSI)
  beta2 <- w

  # ============================================================================
  # STEP 3: Compute restriction matrices K_G1 and K_G2
  # ============================================================================


  # Q_G1 = Gamma1^{-1}A1 C1 (C1'A1 Gamma1^{-1}A1 C1)^{-1} C1'A1
  A1_C1 <- A1 %*% C1
  middle_term_1 <- crossprod(C1, Gamma1_inv_A1) %*% A1_C1

  tryCatch(
    {
      middle_inv_1 <- ginv(middle_term_1)
    },
    error = function(e) {
      stop("Failed to compute restriction matrix for stage 1: ", e$message)
    }
  )

  Q_G1 <- Gamma1_inv_A1 %*% C1 %*% middle_inv_1 %*% crossprod(C1, A1)
  K_G1 <- diag(n1) - Q_G1


  # For stage 2, we need Gamma^{-1}
  Gamma_inv <- matrix(0, nrow = n, ncol = n)
  for (j in seq_len(n)) {
    Gamma_inv[, j] <- cpp_symmetric_solve(Gamma, diag(n)[, j])
  }

  # Q_G2 = Gamma^{-1}Gamma C2 (C2'Gamma Gamma^{-1}Gamma C2)^{-1} C2'Gamma
  # Simplifies to: Gamma^{-1}Gamma C2 (C2'Gamma C2)^{-1} C2'Gamma
  # = C2 (C2'Gamma C2)^{-1} C2'Gamma
  Gamma_C2 <- Gamma %*% C2
  middle_term_2 <- crossprod(C2, Gamma_C2)

  tryCatch(
    {
      middle_inv_2 <- ginv(middle_term_2)
    },
    error = function(e) {
      stop("Failed to compute restriction matrix for stage 2: ", e$message)
    }
  )

  Q_G2 <- C2 %*% middle_inv_2 %*% crossprod(C2, Gamma)
  K_G2 <- diag(n) - Q_G2

  # ============================================================================
  # STEP 4: Apply restrictions
  # ============================================================================

  beta_R1 <- K_G1 %*% beta1
  beta_R2 <- K_G2 %*% beta2

  # ============================================================================
  # STEP 5: Compute phenotypic index for C* adjustment
  # ============================================================================

  P1_inv_G1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    P1_inv_G1[, j] <- cpp_symmetric_solve(P1, G1[, j])
  }
  b1 <- P1_inv_G1 %*% w1

  # Apply same restriction to phenotypic index
  b_R1 <- K_G1 %*% b1

  # ============================================================================
  # STEP 6: Compute correlation and selection intensities
  # ============================================================================

  rho_I1I2 <- .genomic_index_correlation(beta_R1, beta_R2, Gamma1, Gamma, A)

  # Set tau from selection proportion if not provided
  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  if (use_young_method && !is.na(rho_I1I2)) {
    intensities <- tryCatch(
      {
        .young_intensities(selection_proportion, rho_I1I2)
      },
      error = function(e) {
        warning("Young's method failed. Using manual intensities.")
        list(k1 = k1_manual, k2 = k2_manual)
      }
    )
    k1 <- intensities$k1
    k2 <- intensities$k2
  } else {
    k1 <- k1_manual
    k2 <- k2_manual
  }

  # ============================================================================
  # STEP 7: Adjust covariance matrices
  # ============================================================================

  Gamma_star <- .cochran_adjustment_genomic(Gamma, Gamma1, beta_R1, A, k1, tau)

  # Adjust C
  u <- k1 * (k1 - tau)
  b_R1_P1_b_R1 <- cpp_quadratic_form_sym(b_R1, P1)

  if (b_R1_P1_b_R1 > 0) {
    # C* = C - u * C[,1:n1] b_R1 b_R1' C[1:n1,] / (b_R1'P1 b_R1)
    n1 <- nrow(G1)
    C_col_1_n1 <- C[, 1:n1]
    C_row_1_n1 <- C[1:n1, ]

    C_col_b_R1 <- C_col_1_n1 %*% b_R1 # n x 1
    b_R1_C_row <- crossprod(b_R1, C_row_1_n1) # 1 x n
    adjustment_C <- u * (C_col_b_R1 %*% b_R1_C_row) / b_R1_P1_b_R1 # n x n
    C_star <- C - adjustment_C
  } else {
    warning("Invalid phenotypic variance at stage 1. Using unadjusted C.")
    C_star <- C
  }

  # ============================================================================
  # STEP 8: Compute stage metrics
  # ============================================================================

  stage1_metrics <- .genomic_stage_metrics(beta_R1, Gamma1, A1, w1, k1)

  wCw <- cpp_quadratic_form_sym(w, C)
  rho_HI1 <- if (stage1_metrics$beta_Gamma_beta > 0 && wCw > 0) {
    sqrt(stage1_metrics$beta_Gamma_beta / wCw)
  } else {
    NA_real_
  }
  stage1_metrics$rho_HI <- rho_HI1

  stage2_metrics <- .genomic_stage_metrics(beta_R2, Gamma_star, NULL, w, k2)

  wCstarw <- cpp_quadratic_form_sym(w, C_star)
  rho_HI2 <- if (stage2_metrics$beta_Gamma_beta > 0 && wCstarw > 0) {
    sqrt(stage2_metrics$beta_Gamma_beta / wCstarw)
  } else {
    NA_real_
  }
  stage2_metrics$rho_HI <- rho_HI2

  # ============================================================================
  # STEP 9: Create summary
  # ============================================================================

  trait_names_1 <- if (!is.null(colnames(Gamma1))) colnames(Gamma1) else paste0("Trait", 1:n1)
  trait_names_2 <- if (!is.null(colnames(Gamma))) colnames(Gamma) else paste0("Trait", 1:n)

  summary_stage1 <- data.frame(
    Stage = 1,
    Trait = trait_names_1,
    beta_R = round(as.numeric(beta_R1), 4),
    E = round(stage1_metrics$E, 4),
    stringsAsFactors = FALSE
  )

  summary_stage2 <- data.frame(
    Stage = 2,
    Trait = trait_names_2,
    beta_R = round(as.numeric(beta_R2), 4),
    E = round(stage2_metrics$E, 4),
    stringsAsFactors = FALSE
  )

  # ============================================================================
  # STEP 10: Return results
  # ============================================================================

  result <- list(
    beta_R1 = as.numeric(beta_R1),
    beta_R2 = as.numeric(beta_R2),
    beta1 = as.numeric(beta1),
    beta2 = as.numeric(beta2),
    K_G1 = K_G1,
    K_G2 = K_G2,
    stage1_metrics = list(
      R = stage1_metrics$R,
      E = setNames(stage1_metrics$E, trait_names_1),
      rho_HI = stage1_metrics$rho_HI,
      sigma_I = stage1_metrics$sigma_I
    ),
    stage2_metrics = list(
      R = stage2_metrics$R,
      E = setNames(stage2_metrics$E, trait_names_2),
      rho_HI = stage2_metrics$rho_HI,
      sigma_I = stage2_metrics$sigma_I
    ),
    Gamma_star = Gamma_star,
    C_star = C_star,
    rho_I1I2 = rho_I1I2,
    k1 = k1,
    k2 = k2,
    tau = tau,
    selection_proportion = selection_proportion,
    summary_stage1 = summary_stage1,
    summary_stage2 = summary_stage2
  )

  class(result) <- c("mrlgsi", "multistage_genomic_index", "list")

  result
}

# ==============================================================================
# 9.6 MULTISTAGE PREDETERMINED PROPORTIONAL GAIN LGSI (MPPG-LGSI)
# ==============================================================================

#' Multistage Predetermined Proportional Gain Linear Genomic Selection Index (MPPG-LGSI)
#'
#' @description
#' Implements the two-stage Predetermined Proportional Gain LGSI where breeders
#' specify desired proportional gains between traits at each stage using GEBVs.
#'
#' @param Gamma1 GEBV variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param Gamma GEBV variance-covariance matrix for all traits at stage 2 (n x n)
#' @param A1 Covariance matrix between GEBVs and true breeding values for stage 1 (n1 x n1)
#' @param A Covariance matrix between GEBVs and true breeding values for stage 2 (n x n1)
#' @param C Genotypic variance-covariance matrix for all traits (n x n)
#' @param G1 Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param P1 Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param wmat Economic weights vector or matrix (n x k)
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param d1 Vector of desired proportional gains for stage 1 (length n1)
#' @param d2 Vector of desired proportional gains for stage 2 (length n)
#' @param U1 Constraint matrix for stage 1 (n1 x r1), optional
#' @param U2 Constraint matrix for stage 2 (n x r2), optional
#' @param selection_proportion Proportion selected at each stage (default: 0.1)
#' @param use_young_method Logical. Use Young's method for selection intensities (default: FALSE).
#'   Young's method tends to overestimate intensities; manual intensities are recommended.
#' @param k1_manual Manual selection intensity for stage 1
#' @param k2_manual Manual selection intensity for stage 2
#' @param tau Standardized truncation point
#'
#' @return List with components similar to mlgsi, plus:
#'   \itemize{
#'     \item \code{beta_P1} - PPG genomic stage 1 coefficients
#'     \item \code{beta_P2} - PPG genomic stage 2 coefficients
#'     \item \code{b_P1} - PPG phenotypic stage 1 coefficients (used for C* adjustment)
#'     \item \code{theta1} - Proportionality constant for stage 1
#'     \item \code{theta2} - Proportionality constant for stage 2
#'     \item \code{gain_ratios_1} - Achieved gain ratios at stage 1
#'     \item \code{gain_ratios_2} - Achieved gain ratios at stage 2
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The PPG genomic coefficients are:
#' \deqn{\mathbf{\beta}_{P_1} = \mathbf{\beta}_{R_1} + \theta_1 \mathbf{U}_1(\mathbf{U}_1'\mathbf{\Gamma}_1\mathbf{U}_1)^{-1}\mathbf{d}_1}
#' \deqn{\mathbf{\beta}_{P_2} = \mathbf{\beta}_{R_2} + \theta_2 \mathbf{U}_2(\mathbf{U}_2'\mathbf{\Gamma}\mathbf{U}_2)^{-1}\mathbf{d}_2}
#'
#' where proportionality constants are:
#' \deqn{\theta_1 = \frac{\mathbf{d}_1'(\mathbf{U}_1'\mathbf{\Gamma}_1\mathbf{U}_1)^{-1}\mathbf{U}_1'\mathbf{A}_1\mathbf{w}}{\mathbf{d}_1'(\mathbf{U}_1'\mathbf{\Gamma}_1\mathbf{U}_1)^{-1}\mathbf{d}_1}}
#'
#' \strong{Covariance Adjustment:}
#'
#' The genetic covariance matrix \eqn{\mathbf{C}^*} is adjusted using phenotypic PPG coefficients
#' \eqn{\mathbf{b}_{P1} = \mathbf{P}_1^{-1}\mathbf{G}_1\mathbf{P}_1^{-1}\mathbf{d}_1}, which reflect
#' the same proportional gain constraints as the genomic coefficients \eqn{\mathbf{\beta}_{P1}}.
#' This ensures the adjustment reflects the actual PPG selection occurring at stage 1.
#'
#' \strong{Important:} When using custom \code{U1} matrices (subset constraints), the phenotypic
#' proxy \eqn{\mathbf{b}_{P1}} uses the standard Tallis formula (all traits constrained), while
#' the genomic index \eqn{\mathbf{\beta}_{P1}} respects the \code{U1} subset. This may cause
#' \eqn{\mathbf{C}^*} to be slightly over-adjusted. For exact adjustment, use \code{U1 = NULL}
#' (default, all traits constrained). Calculating the exact restricted phenotypic proxy would
#' require implementing the full MPPG-LPSI projection matrix method for the phenotypic coefficients.
#'
#' \strong{Note:} Input covariance matrices (\code{C}, \code{P1}, \code{Gamma1}, \code{Gamma})
#' should be positive definite. Non-positive definite matrices may lead to invalid results or warnings.
#'
#' @references
#' Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing. Chapter 9, Section 9.6.
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-stage proportional gain genomic selection
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' reliability <- 0.7
#' Gamma1 <- reliability * gmat[1:3, 1:3]
#' Gamma <- reliability * gmat
#' A1 <- reliability * gmat[1:3, 1:3]
#' A <- gmat[, 1:3]
#'
#' # Desired proportional gains
#' d1 <- c(2, 1, 1)
#' d2 <- c(3, 2, 1, 1, 1, 0.5, 0.5)
#'
#' weights <- c(10, 8, 6, 4, 3, 2, 1)
#'
#' result <- mppg_lgsi(
#'   Gamma1 = Gamma1, Gamma = Gamma, A1 = A1, A = A,
#'   C = gmat, G1 = gmat[1:3, 1:3], P1 = pmat[1:3, 1:3],
#'   wmat = weights, d1 = d1, d2 = d2
#' )
#' }
mppg_lgsi <- function(Gamma1, Gamma, A1, A, C, G1, P1, wmat, wcol = 1,
                      d1, d2, U1 = NULL, U2 = NULL,
                      selection_proportion = 0.1,
                      use_young_method = FALSE,
                      k1_manual = 2.063,
                      k2_manual = 2.063,
                      tau = NULL) {
  # ============================================================================
  # STEP 1: Input Validation
  # ============================================================================

  Gamma1 <- as.matrix(Gamma1)
  Gamma <- as.matrix(Gamma)
  A1 <- as.matrix(A1)
  d1 <- as.numeric(d1)
  d2 <- as.numeric(d2)

  n1 <- nrow(Gamma1)
  n <- nrow(Gamma)

  if (length(d1) != n1) {
    stop("d1 must have length equal to number of stage 1 traits")
  }

  if (length(d2) != n) {
    stop("d2 must have length equal to number of stage 2 traits")
  }

  # If U matrices not provided, use identity (all traits constrained)
  if (is.null(U1)) {
    U1 <- diag(n1)
  } else {
    U1 <- as.matrix(U1)
    # Check if U1 is not identity matrix
    if (!isTRUE(all.equal(U1, diag(n1), tolerance = 1e-10))) {
      warning("Custom U1 matrix detected. Note: The phenotypic proxy b_P1 uses the ",
        "standard PPG formula (all traits constrained), while beta_P1 respects U1. ",
        "This may cause C* to be slightly over-adjusted. For exact adjustment with ",
        "subset constraints, consider using U1 = Identity (all traits constrained).",
        call. = FALSE
      )
    }
  }

  if (is.null(U2)) {
    U2 <- diag(n)
  } else {
    U2 <- as.matrix(U2)
  }

  # Process weights
  wmat <- as.matrix(wmat)
  if (ncol(wmat) == 1) {
    w <- as.numeric(wmat)
  } else {
    w <- as.numeric(wmat[, wcol])
  }

  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  # ============================================================================
  # STEP 2: Compute restricted coefficients (base for PPG)
  # ============================================================================

  w1 <- w[1:n1]

  # Stage 1: beta1 = Gamma1^{-1} A1 w1
  Gamma1_inv_A1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    Gamma1_inv_A1[, j] <- cpp_symmetric_solve(Gamma1, A1[, j])
  }
  beta1 <- Gamma1_inv_A1 %*% w1

  # For stage 1, compute restriction matrix
  A1_U1 <- A1 %*% U1
  middle_term_1 <- crossprod(U1, Gamma1_inv_A1) %*% A1_U1

  tryCatch(
    {
      middle_inv_1 <- ginv(middle_term_1)
    },
    error = function(e) {
      stop("Failed to compute PPG matrix for stage 1: ", e$message)
    }
  )

  Q_G1 <- Gamma1_inv_A1 %*% U1 %*% middle_inv_1 %*% crossprod(U1, A1)
  K_G1 <- diag(n1) - Q_G1
  beta_R1 <- K_G1 %*% beta1

  # Stage 2: beta2 = w
  beta2 <- w

  # For stage 2, compute restriction matrix
  Gamma_inv <- matrix(0, nrow = n, ncol = n)
  for (j in seq_len(n)) {
    Gamma_inv[, j] <- cpp_symmetric_solve(Gamma, diag(n)[, j])
  }

  Gamma_U2 <- Gamma %*% U2
  middle_term_2 <- crossprod(U2, Gamma_U2)

  tryCatch(
    {
      middle_inv_2 <- ginv(middle_term_2)
    },
    error = function(e) {
      stop("Failed to compute PPG matrix for stage 2: ", e$message)
    }
  )

  Q_G2 <- U2 %*% middle_inv_2 %*% crossprod(U2, Gamma)
  K_G2 <- diag(n) - Q_G2
  beta_R2 <- K_G2 %*% beta2

  # ============================================================================
  # STEP 3: Compute proportionality constants
  # ============================================================================

  # theta1 = d1' (U1'Gamma1 U1)^{-1} U1' A1 w / d1' (U1'Gamma1 U1)^{-1} d1
  U1_Gamma1_U1_inv_d1 <- middle_inv_1 %*% d1
  numerator1 <- as.numeric(crossprod(d1, middle_inv_1) %*% crossprod(U1, A1_U1 %*% w1))
  denominator1 <- as.numeric(crossprod(d1, U1_Gamma1_U1_inv_d1))

  theta1 <- if (abs(denominator1) > 1e-10) numerator1 / denominator1 else 0

  # theta2 = d2' (U2'Gamma U2)^{-1} U2' Gamma w / d2' (U2'Gamma U2)^{-1} d2
  U2_Gamma_U2_inv_d2 <- middle_inv_2 %*% d2
  numerator2 <- as.numeric(crossprod(d2, middle_inv_2) %*% crossprod(U2, Gamma_U2 %*% w))
  denominator2 <- as.numeric(crossprod(d2, U2_Gamma_U2_inv_d2))

  theta2 <- if (abs(denominator2) > 1e-10) numerator2 / denominator2 else 0

  # ============================================================================
  # STEP 4: Compute PPG coefficients
  # ============================================================================

  # beta_P1 = beta_R1 + theta1 * U1 (U1'Gamma1 U1)^{-1} d1
  beta_P1 <- beta_R1 + theta1 * (U1 %*% U1_Gamma1_U1_inv_d1)

  # beta_P2 = beta_R2 + theta2 * U2 (U2'Gamma U2)^{-1} d2
  beta_P2 <- beta_R2 + theta2 * (U2 %*% U2_Gamma_U2_inv_d2)

  # ============================================================================
  # STEP 5: Compute phenotypic PPG coefficients for C* adjustment
  # ============================================================================

  # Compute phenotypic coefficients
  P1_inv_G1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    P1_inv_G1[, j] <- cpp_symmetric_solve(P1, G1[, j])
  }

  # Compute phenotypic PPG coefficients (analogous to beta_P1)
  # b_P1 = P1^{-1} G1 P1^{-1} d1
  P1_inv_d1 <- cpp_symmetric_solve(P1, d1)
  b_P1 <- P1_inv_G1 %*% P1_inv_d1

  # ============================================================================
  # STEP 6: Compute correlation and selection intensities
  # ============================================================================

  rho_I1I2 <- .genomic_index_correlation(beta_P1, beta_P2, Gamma1, Gamma, A)

  # Set tau from selection proportion if not provided
  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  if (use_young_method && !is.na(rho_I1I2)) {
    intensities <- tryCatch(
      {
        .young_intensities(selection_proportion, rho_I1I2)
      },
      error = function(e) {
        warning("Young's method failed. Using manual intensities.")
        list(k1 = k1_manual, k2 = k2_manual)
      }
    )
    k1 <- intensities$k1
    k2 <- intensities$k2
  } else {
    k1 <- k1_manual
    k2 <- k2_manual
  }

  # ============================================================================
  # STEP 7: Adjust covariance matrices
  # ============================================================================

  Gamma_star <- .cochran_adjustment_genomic(Gamma, Gamma1, beta_P1, A, k1, tau)

  # Adjust C using phenotypic PPG coefficients b_P1
  # This ensures C* reflects the actual PPG selection at stage 1
  u <- k1 * (k1 - tau)
  b_P1_P1_b_P1 <- cpp_quadratic_form_sym(b_P1, P1)

  if (b_P1_P1_b_P1 > 0) {
    # C* = C - u * C[,1:n1] b_P1 b_P1' C[1:n1,] / (b_P1'P1 b_P1)
    n1 <- nrow(G1)
    C_col_1_n1 <- C[, 1:n1]
    C_row_1_n1 <- C[1:n1, ]

    C_col_b_P1 <- C_col_1_n1 %*% b_P1 # n x 1
    b_P1_C_row <- crossprod(b_P1, C_row_1_n1) # 1 x n
    adjustment_C <- u * (C_col_b_P1 %*% b_P1_C_row) / b_P1_P1_b_P1 # n x n
    C_star <- C - adjustment_C

    # Check if C_star is still positive definite
    min_eig_C_star <- tryCatch(
      {
        min(eigen(C_star, symmetric = TRUE, only.values = TRUE)$values)
      },
      error = function(e) -Inf
    )

    if (min_eig_C_star <= 1e-8) {
      warning(
        "C* is not positive definite (min eigenvalue = ", round(min_eig_C_star, 6),
        "). This may indicate issues with input covariance matrices. ",
        "Consider checking that C and P1 are positive definite."
      )
    }
  } else {
    warning("Invalid phenotypic variance at stage 1. Using unadjusted C.")
    C_star <- C
  }

  # ============================================================================
  # STEP 8: Compute stage metrics
  # ============================================================================

  stage1_metrics <- .genomic_stage_metrics(beta_P1, Gamma1, A1, w1, k1)

  wCw <- cpp_quadratic_form_sym(w, C)
  rho_HI1 <- if (stage1_metrics$beta_Gamma_beta > 0 && wCw > 0) {
    sqrt(stage1_metrics$beta_Gamma_beta / wCw)
  } else {
    NA_real_
  }
  stage1_metrics$rho_HI <- rho_HI1

  stage2_metrics <- .genomic_stage_metrics(beta_P2, Gamma_star, NULL, w, k2)

  wCstarw <- cpp_quadratic_form_sym(w, C_star)
  rho_HI2 <- if (stage2_metrics$beta_Gamma_beta > 0 && wCstarw > 0) {
    sqrt(stage2_metrics$beta_Gamma_beta / wCstarw)
  } else {
    NA_real_
  }
  stage2_metrics$rho_HI <- rho_HI2

  # ============================================================================
  # STEP 9: Compute gain ratios
  # ============================================================================

  gain_ratios_1 <- stage1_metrics$E / d1
  gain_ratios_1[!is.finite(gain_ratios_1)] <- NA_real_

  gain_ratios_2 <- stage2_metrics$E / d2
  gain_ratios_2[!is.finite(gain_ratios_2)] <- NA_real_

  # ============================================================================
  # STEP 10: Create summary
  # ============================================================================

  trait_names_1 <- if (!is.null(colnames(Gamma1))) colnames(Gamma1) else paste0("Trait", 1:n1)
  trait_names_2 <- if (!is.null(colnames(Gamma))) colnames(Gamma) else paste0("Trait", 1:n)

  summary_stage1 <- data.frame(
    Stage = 1,
    Trait = trait_names_1,
    beta_P = round(as.numeric(beta_P1), 4),
    d = round(d1, 4),
    E = round(stage1_metrics$E, 4),
    Ratio = round(gain_ratios_1, 4),
    stringsAsFactors = FALSE
  )

  summary_stage2 <- data.frame(
    Stage = 2,
    Trait = trait_names_2,
    beta_P = round(as.numeric(beta_P2), 4),
    d = round(d2, 4),
    E = round(stage2_metrics$E, 4),
    Ratio = round(gain_ratios_2, 4),
    stringsAsFactors = FALSE
  )

  # ============================================================================
  # STEP 11: Return results
  # ============================================================================

  result <- list(
    beta_P1 = as.numeric(beta_P1),
    beta_P2 = as.numeric(beta_P2),
    beta_R1 = as.numeric(beta_R1),
    beta_R2 = as.numeric(beta_R2),
    b_P1 = as.numeric(b_P1), # Phenotypic PPG coefficients used for C* adjustment
    theta1 = theta1,
    theta2 = theta2,
    d1 = d1,
    d2 = d2,
    gain_ratios_1 = setNames(gain_ratios_1, trait_names_1),
    gain_ratios_2 = setNames(gain_ratios_2, trait_names_2),
    stage1_metrics = list(
      R = stage1_metrics$R,
      E = setNames(stage1_metrics$E, trait_names_1),
      rho_HI = stage1_metrics$rho_HI,
      sigma_I = stage1_metrics$sigma_I
    ),
    stage2_metrics = list(
      R = stage2_metrics$R,
      E = setNames(stage2_metrics$E, trait_names_2),
      rho_HI = stage2_metrics$rho_HI,
      sigma_I = stage2_metrics$sigma_I
    ),
    Gamma_star = Gamma_star,
    C_star = C_star,
    rho_I1I2 = rho_I1I2,
    k1 = k1,
    k2 = k2,
    tau = tau,
    selection_proportion = selection_proportion,
    summary_stage1 = summary_stage1,
    summary_stage2 = summary_stage2
  )

  class(result) <- c("mppg_lgsi", "multistage_genomic_index", "list")

  result
}
