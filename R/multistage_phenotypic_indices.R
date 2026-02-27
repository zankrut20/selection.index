#' Multistage Linear Phenotypic Selection Indices (Chapter 9)
#' @name multistage_phenotypic_indices
#'
#' @description
#' Implements multistage phenotypic selection index methods from Chapter 9.
#' These methods involve selection across multiple stages with covariance
#' adjustments due to selection effects using Cochran/Cunningham's method.
#'
#' Methods included:
#' - MLPSI: Multistage Linear Phenotypic Selection Index (Section 9.1)
#' - MRLPSI: Multistage Restricted Linear Phenotypic Selection Index (Section 9.2)
#' - MPPG-LPSI: Multistage Predetermined Proportional Gain LPSI (Section 9.3)
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
#' @importFrom stats setNames pnorm qnorm dnorm
#' @importFrom MASS ginv
NULL


#' Compute Cochran/Cunningham covariance adjustment
#' @keywords internal
#' @noRd
.cochran_adjustment <- function(cov_xy, cov_xw, cov_yw, var_w, k1, tau) {
  u <- k1 * (k1 - tau)
  cov_xy - u * (cov_xw * cov_yw) / var_w
}

#' Adjust phenotypic covariance matrix after stage 1 selection
#' @param stage1_indices Indices of stage 1 traits in the full matrix (default: 1:nrow(P1))
#' @keywords internal
#' @noRd
.adjust_phenotypic_matrix <- function(P, P1, b1, k1, tau, stage1_indices = NULL) {
  u <- k1 * (k1 - tau)
  b1Pb1 <- cpp_quadratic_form_sym(b1, P1)

  if (b1Pb1 <= 0) {
    warning("Invalid variance at stage 1 (b1'P1b1 <= 0). Returning unadjusted matrix.")
    return(P)
  }

  n1 <- nrow(P1)
  if (is.null(stage1_indices)) {
    stage1_indices <- 1:n1
  }

  if (length(stage1_indices) != n1) {
    stop("Length of stage1_indices must match nrow(P1)")
  }

  P_col_stage1 <- P[, stage1_indices]
  P_row_stage1 <- P[stage1_indices, ]

  P_col_b1 <- P_col_stage1 %*% b1 # n x 1
  b1_P_row <- crossprod(b1, P_row_stage1) # 1 x n

  outer_prod <- P_col_b1 %*% b1_P_row

  scalar_denom <- as.numeric(b1Pb1)

  temp1 <- u * outer_prod
  adjustment <- temp1 / scalar_denom

  P_star <- P - adjustment
  P_star <- as.matrix(P_star)

  P_star
}

#' Adjust genotypic covariance matrix after stage 1 selection
#' @param stage1_indices Indices of stage 1 traits in the full matrix (default: 1:nrow(G1))
#' @keywords internal
#' @noRd
.adjust_genotypic_matrix <- function(C, G1, b1, k1, tau, P1, stage1_indices = NULL) {
  u <- k1 * (k1 - tau)
  b1Pb1 <- cpp_quadratic_form_sym(b1, P1)

  if (b1Pb1 <= 0) {
    warning("Invalid variance at stage 1 (b1'P1b1 <= 0). Returning unadjusted matrix.")
    return(C)
  }

  n1 <- nrow(G1)
  if (is.null(stage1_indices)) {
    stage1_indices <- 1:n1
  }

  if (length(stage1_indices) != n1) {
    stop("Length of stage1_indices must match nrow(G1)")
  }



  C_col_stage1 <- C[, stage1_indices]
  C_row_stage1 <- C[stage1_indices, ]

  C_stage1_block <- C[stage1_indices, stage1_indices]
  if (!isTRUE(all.equal(G1, C_stage1_block, tolerance = 1e-6))) {
    warning(
      "G1 does not match C[stage1_indices, stage1_indices]. ",
      "This may indicate incorrect stage1_indices or misaligned matrices."
    )
  }

  C_col_b1 <- C_col_stage1 %*% b1 # n x 1
  b1_C_row <- crossprod(b1, C_row_stage1) # 1 x n
  adjustment <- u * (C_col_b1 %*% b1_C_row) / b1Pb1 # n x n
  C_star <- C - adjustment

  C_star
}

#' Compute correlation between indices at two stages
#' @param stage1_indices Indices of stage 1 traits in the full matrix (default: 1:nrow(P1))
#' @keywords internal
#' @noRd
.index_correlation <- function(b1, b2, P1, P, stage1_indices = NULL) {
  n1 <- nrow(P1)
  if (is.null(stage1_indices)) {
    stage1_indices <- 1:n1
  }

  b1Pb1 <- cpp_quadratic_form_sym(b1, P1)
  b2Pb2 <- cpp_quadratic_form_sym(b2, P)

  if (b1Pb1 <= 0 || b2Pb2 <= 0) {
    warning("Invalid variance for correlation calculation.")
    return(NA_real_)
  }

  P_stage1_all <- P[stage1_indices, ]

  numerator <- as.numeric(crossprod(b1, P_stage1_all %*% b2))
  denominator <- sqrt(b1Pb1) * sqrt(b2Pb2)

  rho_12 <- numerator / denominator
  rho_12
}

#' Young's method for selection intensities
#' @keywords internal
#' @noRd
.young_selection_intensities <- function(p, rho_12) {

  warning("Young's method tends to overestimate selection intensities. ",
    "Consider using manual intensities (use_young_method = FALSE) for more conservative estimates.",
    call. = FALSE
  )

  if (p <= 0 || p >= 1) {
    stop("Selection proportion p must be between 0 and 1")
  }

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

#' Compute stage metrics for multistage indices
#' @keywords internal
#' @noRd
.stage_metrics <- function(b, P, G, w, k) {
  b <- as.numeric(b)
  bPb <- cpp_quadratic_form_sym(b, P)
  bGb <- cpp_quadratic_form_sym(b, G)

  sigma_I <- if (bPb > 0) sqrt(bPb) else NA_real_

  R <- if (!is.na(sigma_I)) k * sigma_I else NA_real_

  Gb <- G %*% b
  E <- if (!is.na(sigma_I) && sigma_I > 0) k * Gb / sigma_I else rep(NA_real_, length(Gb))

  wCw <- cpp_quadratic_form_sym(w, G)
  rho_H <- if (!is.na(bPb) && bPb > 0 && wCw > 0) sqrt(bPb / wCw) else NA_real_

  list(
    sigma_I = sigma_I,
    R = R,
    E = as.numeric(E),
    rho_H = rho_H,
    bPb = bPb,
    bGb = bGb
  )
}


#' Multistage Linear Phenotypic Selection Index (MLPSI)
#'
#' @description
#' Implements the two-stage Linear Phenotypic Selection Index where selection
#' occurs at two different stages with covariance adjustments due to selection
#' effects using Cochran/Cunningham's method.
#'
#' @param P1 Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param P Phenotypic variance-covariance matrix for all traits at stage 2 (n x n)
#' @param G1 Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param C Genotypic variance-covariance matrix for all traits (n x n)
#' @param wmat Economic weights vector or matrix (n x k)
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param stage1_indices Integer vector specifying which traits (columns of P and C) correspond
#'   to stage 1. Default is 1:nrow(P1), assuming stage 1 traits are the first n1 traits.
#'   Use this to specify non-contiguous traits, e.g., c(1, 3, 5) for traits 1, 3, and 5.
#' @param selection_proportion Proportion selected at each stage (default: 0.1)
#' @param use_young_method Logical. Use Young's method for selection intensities (default: FALSE).
#'   Young's method tends to overestimate intensities; manual intensities are recommended.
#' @param k1_manual Manual selection intensity for stage 1 (used if use_young_method = FALSE)
#' @param k2_manual Manual selection intensity for stage 2 (used if use_young_method = FALSE)
#' @param tau Standardized truncation point (default: computed from selection_proportion)
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{b1} - Stage 1 index coefficients
#'     \item \code{b2} - Stage 2 index coefficients
#'     \item \code{stage1_metrics} - List with stage 1 metrics (R1, E1, rho_H1)
#'     \item \code{stage2_metrics} - List with stage 2 metrics (R2, E2, rho_H2)
#'     \item \code{P_star} - Adjusted phenotypic covariance matrix at stage 2
#'     \item \code{C_star} - Adjusted genotypic covariance matrix at stage 2
#'     \item \code{rho_12} - Correlation between stage 1 and stage 2 indices
#'     \item \code{k1} - Selection intensity at stage 1
#'     \item \code{k2} - Selection intensity at stage 2
#'     \item \code{summary_stage1} - Data frame with stage 1 summary
#'     \item \code{summary_stage2} - Data frame with stage 2 summary
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' Stage 1 index coefficients:
#' \deqn{\mathbf{b}_1 = \mathbf{P}_1^{-1}\mathbf{G}_1\mathbf{w}}
#'
#' Stage 2 index coefficients:
#' \deqn{\mathbf{b}_2 = \mathbf{P}^{-1}\mathbf{G}\mathbf{w}}
#'
#' Adjusted phenotypic covariance matrix (Cochran/Cunningham):
#' \deqn{\mathbf{P}^* = \mathbf{P} - u \frac{Cov(\mathbf{y},\mathbf{x}_1)\mathbf{b}_1\mathbf{b}_1'Cov(\mathbf{x}_1,\mathbf{y})}{\mathbf{b}_1'\mathbf{P}_1\mathbf{b}_1}}
#'
#' Adjusted genotypic covariance matrix:
#' \deqn{\mathbf{C}^* = \mathbf{C} - u \frac{\mathbf{G}_1'\mathbf{b}_1\mathbf{b}_1'\mathbf{G}_1}{\mathbf{b}_1'\mathbf{P}_1\mathbf{b}_1}}
#'
#' where \eqn{u = k_1(k_1 - \tau)}
#'
#' Selection response: \eqn{R_1 = k_1 \sqrt{\mathbf{b}_1'\mathbf{P}_1\mathbf{b}_1}},
#' \eqn{R_2 = k_2 \sqrt{\mathbf{b}_2'\mathbf{P}^*\mathbf{b}_2}}
#'
#' @references
#' Cochran, W. G. (1951). Improvement by means of selection.
#' Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability, 449-470.
#'
#' Cunningham, E. P. (1975). Multi-stage index selection.
#' Theoretical and Applied Genetics, 46(2), 55-61.
#'
#' Young, S. S. Y. (1964). Multi-stage selection for genetic gain.
#' Heredity, 19(1), 131-144.
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-stage selection example
#' # Stage 1: Select based on 3 traits
#' # Stage 2: Select based on all 7 traits
#'
#' # Compute variance-covariance matrices
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Stage 1 uses first 3 traits
#' P1 <- pmat[1:3, 1:3]
#' G1 <- gmat[1:3, 1:3]
#'
#' # Stage 2 uses all 7 traits
#' P <- pmat
#' C <- gmat
#'
#' # Economic weights
#' weights <- c(10, 8, 6, 4, 3, 2, 1)
#'
#' # Run MLPSI (default: stage1_indices = 1:3)
#' result <- mlpsi(
#'   P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
#'   selection_proportion = 0.1
#' )
#'
#' # Or with non-contiguous traits (e.g., traits 1, 3, 5 at stage 1):
#' # P1 <- pmat[c(1,3,5), c(1,3,5)]
#' # G1 <- gmat[c(1,3,5), c(1,3,5)]
#' # result <- mlpsi(P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
#' #                 stage1_indices = c(1, 3, 5))
#'
#' print(result$summary_stage1)
#' print(result$summary_stage2)
#' }
mlpsi <- function(P1, P, G1, C, wmat, wcol = 1,
                  stage1_indices = NULL,
                  selection_proportion = 0.1,
                  use_young_method = FALSE,
                  k1_manual = 2.063,
                  k2_manual = 2.063,
                  tau = NULL) {

  P1 <- as.matrix(P1)
  P <- as.matrix(P)
  G1 <- as.matrix(G1)
  C <- as.matrix(C)

  n1 <- nrow(P1)
  n <- nrow(P)

  if (is.null(stage1_indices)) {
    stage1_indices <- 1:n1
  }

  if (length(stage1_indices) != n1) {
    stop("Length of stage1_indices must equal nrow(P1)")
  }
  if (any(stage1_indices < 1) || any(stage1_indices > n)) {
    stop("stage1_indices must be between 1 and nrow(P)")
  }

  wmat <- as.matrix(wmat)
  if (ncol(wmat) == 1) {
    w <- as.numeric(wmat)
  } else {
    w <- as.numeric(wmat[, wcol])
  }

  if (length(w) != n) {
    stop("Weight vector length must match number of traits at stage 2")
  }

  w1 <- w[stage1_indices]
  P1_inv_G1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    P1_inv_G1[, j] <- cpp_symmetric_solve(P1, G1[, j])
  }
  b1 <- P1_inv_G1 %*% w1


  P_inv_G <- matrix(0, nrow = n, ncol = n)
  for (j in seq_len(n)) {
    P_inv_G[, j] <- cpp_symmetric_solve(P, C[, j])
  }
  b2 <- P_inv_G %*% w


  rho_12 <- .index_correlation(b1, b2, P1, P, stage1_indices)


  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  if (use_young_method && !is.na(rho_12)) {
    intensities <- tryCatch(
      {
        .young_selection_intensities(selection_proportion, rho_12)
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


  P_star <- .adjust_phenotypic_matrix(P, P1, b1, k1, tau, stage1_indices)
  C_star <- .adjust_genotypic_matrix(C, G1, b1, k1, tau, P1, stage1_indices)


  stage1_metrics <- .stage_metrics(b1, P1, G1, w1, k1)
  stage2_metrics <- .stage_metrics(b2, P_star, C_star, w, k2)


  trait_names_1 <- if (!is.null(colnames(P1))) colnames(P1) else paste0("Trait", 1:n1)
  trait_names_2 <- if (!is.null(colnames(P))) colnames(P) else paste0("Trait", 1:n)

  summary_stage1 <- data.frame(
    Stage = 1,
    Trait = trait_names_1,
    b = round(as.numeric(b1), 4),
    E = round(stage1_metrics$E, 4),
    stringsAsFactors = FALSE
  )

  summary_stage2 <- data.frame(
    Stage = 2,
    Trait = trait_names_2,
    b = round(as.numeric(b2), 4),
    E = round(stage2_metrics$E, 4),
    stringsAsFactors = FALSE
  )


  result <- list(
    b1 = as.numeric(b1),
    b2 = as.numeric(b2),
    stage1_metrics = list(
      R = stage1_metrics$R,
      E = setNames(stage1_metrics$E, trait_names_1),
      rho_H = stage1_metrics$rho_H,
      sigma_I = stage1_metrics$sigma_I
    ),
    stage2_metrics = list(
      R = stage2_metrics$R,
      E = setNames(stage2_metrics$E, trait_names_2),
      rho_H = stage2_metrics$rho_H,
      sigma_I = stage2_metrics$sigma_I
    ),
    P_star = P_star,
    C_star = C_star,
    rho_12 = rho_12,
    k1 = k1,
    k2 = k2,
    tau = tau,
    selection_proportion = selection_proportion,
    summary_stage1 = summary_stage1,
    summary_stage2 = summary_stage2
  )

  class(result) <- c("mlpsi", "multistage_index", "list")

  result
}


#' Multistage Restricted Linear Phenotypic Selection Index (MRLPSI)
#'
#' @description
#' Implements the two-stage Restricted Linear Phenotypic Selection Index where
#' certain traits are constrained to have zero genetic gain at each stage.
#'
#' @param P1 Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param P Phenotypic variance-covariance matrix for all traits at stage 2 (n x n)
#' @param G1 Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param C Genotypic variance-covariance matrix for all traits (n x n)
#' @param wmat Economic weights vector or matrix (n x k)
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param C1 Constraint matrix for stage 1 (n1 x r1)
#' @param C2 Constraint matrix for stage 2 (n x r2)
#' @param stage1_indices Integer vector specifying which traits correspond to stage 1 (default: 1:nrow(P1))
#' @param selection_proportion Proportion selected at each stage (default: 0.1)
#' @param use_young_method Logical. Use Young's method for selection intensities (default: FALSE).
#'   Young's method tends to overestimate intensities; manual intensities are recommended.
#' @param k1_manual Manual selection intensity for stage 1
#' @param k2_manual Manual selection intensity for stage 2
#' @param tau Standardized truncation point
#'
#' @return List with components similar to mlpsi, plus:
#'   \itemize{
#'     \item \code{b_R1} - Restricted stage 1 coefficients
#'     \item \code{b_R2} - Restricted stage 2 coefficients
#'     \item \code{K1} - Restriction matrix for stage 1
#'     \item \code{K2} - Restriction matrix for stage 2
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The restricted coefficients are computed as:
#' \deqn{\mathbf{b}_{R_1} = \mathbf{K}_1 \mathbf{b}_1}
#' \deqn{\mathbf{b}_{R_2} = \mathbf{K}_2 \mathbf{b}_2}
#'
#' where \eqn{\mathbf{K}_1 = \mathbf{I}_1 - \mathbf{Q}_1} and \eqn{\mathbf{K}_2 = \mathbf{I}_2 - \mathbf{Q}_2}
#'
#' and \eqn{\mathbf{Q}_i = \mathbf{P}_i^{-1}\mathbf{G}_i\mathbf{C}_i(\mathbf{C}_i'\mathbf{G}_i\mathbf{P}_i^{-1}\mathbf{G}_i\mathbf{C}_i)^{-1}\mathbf{C}_i'\mathbf{G}_i}
#'
#' @references
#' Kempthorne, O., & Nordskog, A. W. (1959). Restricted selection indices.
#' Biometrics, 15(1), 10-19.
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-stage restricted selection
#' # Restrict trait 1 at stage 1, traits 1 and 3 at stage 2
#'
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' P1 <- pmat[1:3, 1:3]
#' G1 <- gmat[1:3, 1:3]
#' P <- pmat
#' C <- gmat
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
#' result <- mrlpsi(
#'   P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
#'   C1 = C1, C2 = C2
#' )
#' }
mrlpsi <- function(P1, P, G1, C, wmat, wcol = 1,
                   C1, C2,
                   stage1_indices = NULL,
                   selection_proportion = 0.1,
                   use_young_method = FALSE,
                   k1_manual = 2.063,
                   k2_manual = 2.063,
                   tau = NULL) {

  P1 <- as.matrix(P1)
  P <- as.matrix(P)
  G1 <- as.matrix(G1)
  C <- as.matrix(C)
  C1 <- as.matrix(C1)
  C2 <- as.matrix(C2)

  n1 <- nrow(P1)
  n <- nrow(P)

  if (is.null(stage1_indices)) {
    stage1_indices <- 1:n1
  }

  if (length(stage1_indices) != n1) {
    stop("Length of stage1_indices must equal nrow(P1)")
  }

  wmat <- as.matrix(wmat)
  if (ncol(wmat) == 1) {
    w <- as.numeric(wmat)
  } else {
    w <- as.numeric(wmat[, wcol])
  }

  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }

  w1 <- w[stage1_indices]


  P1_inv_G1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    P1_inv_G1[, j] <- cpp_symmetric_solve(P1, G1[, j])
  }
  b1 <- P1_inv_G1 %*% w1

  P_inv_C <- matrix(0, nrow = n, ncol = n)
  for (j in seq_len(n)) {
    P_inv_C[, j] <- cpp_symmetric_solve(P, C[, j])
  }
  b2 <- P_inv_C %*% w


  G1_C1 <- G1 %*% C1
  middle_term_1 <- crossprod(C1, P1_inv_G1) %*% G1_C1

  tryCatch(
    {
      middle_inv_1 <- ginv(middle_term_1)
    },
    error = function(e) {
      stop("Failed to compute restriction matrix for stage 1: ", e$message)
    }
  )

  Q1 <- P1_inv_G1 %*% C1 %*% middle_inv_1 %*% crossprod(C1, G1)
  K1 <- diag(n1) - Q1


  C_C2 <- C %*% C2
  middle_term_2 <- crossprod(C2, P_inv_C) %*% C_C2

  tryCatch(
    {
      middle_inv_2 <- ginv(middle_term_2)
    },
    error = function(e) {
      stop("Failed to compute restriction matrix for stage 2: ", e$message)
    }
  )

  Q2 <- P_inv_C %*% C2 %*% middle_inv_2 %*% crossprod(C2, C)
  K2 <- diag(n) - Q2


  b_R1 <- K1 %*% b1
  b_R2 <- K2 %*% b2


  rho_12 <- .index_correlation(b_R1, b_R2, P1, P, stage1_indices)

  if (use_young_method && !is.na(rho_12)) {
    intensities <- tryCatch(
      {
        .young_selection_intensities(selection_proportion, rho_12)
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


  P_star <- .adjust_phenotypic_matrix(P, P1, b_R1, k1, tau, stage1_indices)
  C_star <- .adjust_genotypic_matrix(C, G1, b_R1, k1, tau, P1, stage1_indices)


  stage1_metrics <- .stage_metrics(b_R1, P1, G1, w1, k1)
  stage2_metrics <- .stage_metrics(b_R2, P_star, C_star, w, k2)


  trait_names_1 <- if (!is.null(colnames(P1))) colnames(P1) else paste0("Trait", 1:n1)
  trait_names_2 <- if (!is.null(colnames(P))) colnames(P) else paste0("Trait", 1:n)

  summary_stage1 <- data.frame(
    Stage = 1,
    Trait = trait_names_1,
    b_R = round(as.numeric(b_R1), 4),
    E = round(stage1_metrics$E, 4),
    stringsAsFactors = FALSE
  )

  summary_stage2 <- data.frame(
    Stage = 2,
    Trait = trait_names_2,
    b_R = round(as.numeric(b_R2), 4),
    E = round(stage2_metrics$E, 4),
    stringsAsFactors = FALSE
  )


  result <- list(
    b_R1 = as.numeric(b_R1),
    b_R2 = as.numeric(b_R2),
    b1 = as.numeric(b1),
    b2 = as.numeric(b2),
    K1 = K1,
    K2 = K2,
    stage1_metrics = list(
      R = stage1_metrics$R,
      E = setNames(stage1_metrics$E, trait_names_1),
      rho_H = stage1_metrics$rho_H,
      sigma_I = stage1_metrics$sigma_I
    ),
    stage2_metrics = list(
      R = stage2_metrics$R,
      E = setNames(stage2_metrics$E, trait_names_2),
      rho_H = stage2_metrics$rho_H,
      sigma_I = stage2_metrics$sigma_I
    ),
    P_star = P_star,
    C_star = C_star,
    rho_12 = rho_12,
    k1 = k1,
    k2 = k2,
    tau = tau,
    selection_proportion = selection_proportion,
    summary_stage1 = summary_stage1,
    summary_stage2 = summary_stage2
  )

  class(result) <- c("mrlpsi", "multistage_index", "list")

  result
}


#' Multistage Predetermined Proportional Gain Linear Phenotypic Selection Index (MPPG-LPSI)
#'
#' @description
#' Implements the two-stage Predetermined Proportional Gain LPSI where breeders
#' specify desired proportional gains between traits at each stage.
#'
#' @param P1 Phenotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param P Phenotypic variance-covariance matrix for all traits at stage 2 (n x n)
#' @param G1 Genotypic variance-covariance matrix for stage 1 traits (n1 x n1)
#' @param C Genotypic variance-covariance matrix for all traits (n x n)
#' @param wmat Economic weights vector or matrix (n x k)
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param d1 Vector of desired proportional gains for stage 1 (length n1)
#' @param d2 Vector of desired proportional gains for stage 2 (length n)
#' @param stage1_indices Integer vector specifying which traits correspond to stage 1 (default: 1:nrow(P1))
#' @param selection_proportion Proportion selected at each stage (default: 0.1)
#' @param use_young_method Logical. Use Young's method for selection intensities (default: FALSE).
#'   Young's method tends to overestimate intensities; manual intensities are recommended.
#' @param k1_manual Manual selection intensity for stage 1
#' @param k2_manual Manual selection intensity for stage 2
#' @param tau Standardized truncation point
#'
#' @return List with components similar to mlpsi, plus:
#'   \itemize{
#'     \item \code{b_M1} - PPG stage 1 coefficients
#'     \item \code{b_M2} - PPG stage 2 coefficients
#'     \item \code{b_R1} - Restricted stage 1 coefficients
#'     \item \code{b_R2} - Restricted stage 2 coefficients
#'     \item \code{K_M1} - PPG projection matrix for stage 1
#'     \item \code{K_M2} - PPG projection matrix for stage 2
#'     \item \code{theta1} - Proportionality constant for stage 1
#'     \item \code{theta2} - Proportionality constant for stage 2
#'     \item \code{gain_ratios_1} - Achieved gain ratios at stage 1
#'     \item \code{gain_ratios_2} - Achieved gain ratios at stage 2
#'   }
#'
#' @details
#' \strong{Mathematical Formulation (Chapter 9.3.1, Eq 9.17):}
#'
#' The PPG coefficients are computed using the projection matrix method:
#' \deqn{\mathbf{b}_{M_1} = \mathbf{b}_{R_1} + \theta_1 \mathbf{U}_1(\mathbf{U}_1'\mathbf{G}_1\mathbf{P}_1^{-1}\mathbf{G}_1\mathbf{U}_1)^{-1}\mathbf{d}_1}
#' \deqn{\mathbf{b}_{M_2} = \mathbf{b}_{R_2} + \theta_2 \mathbf{U}_2(\mathbf{U}_2'\mathbf{C}\mathbf{P}^{-1}\mathbf{C}\mathbf{U}_2)^{-1}\mathbf{d}_2}
#'
#' where:
#' \itemize{
#'   \item \eqn{\mathbf{b}_{R_i} = \mathbf{K}_{M_i}\mathbf{b}_i} are restricted coefficients
#'   \item \eqn{\mathbf{K}_{M_i} = \mathbf{I} - \mathbf{Q}_{M_i}} is the projection matrix
#'   \item \eqn{\theta_i} is the proportionality constant computed from \eqn{\mathbf{d}_i}
#'   \item \eqn{\mathbf{U}_i = \mathbf{I}} (all traits constrained)
#' }
#'
#' \deqn{\mathbf{b}_{M_1} = \mathbf{K}_{M_1} \mathbf{b}_1}
#' \deqn{\mathbf{b}_{M_2} = \mathbf{K}_{M_2} \mathbf{b}_2}
#'
#' where \eqn{\mathbf{K}_{M_i}} is computed to achieve proportional gains specified by \eqn{\mathbf{d}_i}
#'
#' @references
#' Tallis, G. M. (1962). A selection index for optimum genotype.
#' Biometrics, 18(1), 120-122.
#'
#' @export
#' @examples
#' \dontrun{
#' # Two-stage proportional gain selection
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' P1 <- pmat[1:3, 1:3]
#' G1 <- gmat[1:3, 1:3]
#' P <- pmat
#' C <- gmat
#'
#' # Desired proportional gains
#' d1 <- c(2, 1, 1) # Trait 1 gains twice as much at stage 1
#' d2 <- c(3, 2, 1, 1, 1, 0.5, 0.5) # Different proportions at stage 2
#'
#' weights <- c(10, 8, 6, 4, 3, 2, 1)
#'
#' result <- mppg_lpsi(
#'   P1 = P1, P = P, G1 = G1, C = C, wmat = weights,
#'   d1 = d1, d2 = d2
#' )
#' }
mppg_lpsi <- function(P1, P, G1, C, wmat, wcol = 1,
                      d1, d2,
                      stage1_indices = NULL,
                      selection_proportion = 0.1,
                      use_young_method = FALSE,
                      k1_manual = 2.063,
                      k2_manual = 2.063,
                      tau = NULL) {

  P1 <- as.matrix(P1)
  P <- as.matrix(P)
  G1 <- as.matrix(G1)
  C <- as.matrix(C)
  d1 <- as.numeric(d1)
  d2 <- as.numeric(d2)

  n1 <- nrow(P1)
  n <- nrow(P)

  if (is.null(stage1_indices)) {
    stage1_indices <- 1:n1
  }

  if (length(stage1_indices) != n1) {
    stop("Length of stage1_indices must equal nrow(P1)")
  }

  if (length(d1) != n1) {
    stop("d1 must have length equal to number of stage 1 traits")
  }

  if (length(d2) != n) {
    stop("d2 must have length equal to number of stage 2 traits")
  }

  wmat <- as.matrix(wmat)
  if (ncol(wmat) == 1) {
    w <- as.numeric(wmat)
  } else {
    w <- as.numeric(wmat[, wcol])
  }

  if (is.null(tau)) {
    tau <- qnorm(1 - selection_proportion)
  }


  w1 <- w[stage1_indices]

  P1_inv_G1 <- matrix(0, nrow = n1, ncol = n1)
  for (j in seq_len(n1)) {
    P1_inv_G1[, j] <- cpp_symmetric_solve(P1, G1[, j])
  }
  b1 <- P1_inv_G1 %*% w1

  P_inv_C <- matrix(0, nrow = n, ncol = n)
  for (j in seq_len(n)) {
    P_inv_C[, j] <- cpp_symmetric_solve(P, C[, j])
  }
  b2 <- P_inv_C %*% w


  U1 <- diag(n1)

  G1_U1 <- G1 %*% U1
  middle_term_1 <- crossprod(U1, P1_inv_G1) %*% G1_U1

  tryCatch(
    {
      middle_inv_1 <- ginv(middle_term_1)
    },
    error = function(e) {
      stop("Failed to compute PPG matrix for stage 1: ", e$message)
    }
  )

  Q_M1 <- P1_inv_G1 %*% U1 %*% middle_inv_1 %*% crossprod(U1, G1)
  K_M1 <- diag(n1) - Q_M1
  b_R1 <- K_M1 %*% b1

  U1_G1_P1inv_G1_U1_inv_d1 <- middle_inv_1 %*% d1
  numerator1 <- as.numeric(crossprod(d1, middle_inv_1) %*% crossprod(U1, G1_U1 %*% w1))
  denominator1 <- as.numeric(crossprod(d1, U1_G1_P1inv_G1_U1_inv_d1))

  theta1 <- if (abs(denominator1) > 1e-10) numerator1 / denominator1 else 0

  b_M1 <- b_R1 + theta1 * (U1 %*% U1_G1_P1inv_G1_U1_inv_d1)

  U2 <- diag(n)

  C_U2 <- C %*% U2
  middle_term_2 <- crossprod(U2, P_inv_C) %*% C_U2

  tryCatch(
    {
      middle_inv_2 <- ginv(middle_term_2)
    },
    error = function(e) {
      stop("Failed to compute PPG matrix for stage 2: ", e$message)
    }
  )

  Q_M2 <- P_inv_C %*% U2 %*% middle_inv_2 %*% crossprod(U2, C)
  K_M2 <- diag(n) - Q_M2
  b_R2 <- K_M2 %*% b2

  U2_C_Pinv_C_U2_inv_d2 <- middle_inv_2 %*% d2
  numerator2 <- as.numeric(crossprod(d2, middle_inv_2) %*% crossprod(U2, C_U2 %*% w))
  denominator2 <- as.numeric(crossprod(d2, U2_C_Pinv_C_U2_inv_d2))

  theta2 <- if (abs(denominator2) > 1e-10) numerator2 / denominator2 else 0

  b_M2 <- b_R2 + theta2 * (U2 %*% U2_C_Pinv_C_U2_inv_d2)


  rho_12 <- .index_correlation(b_M1, b_M2, P1, P, stage1_indices)

  if (use_young_method && !is.na(rho_12)) {
    intensities <- tryCatch(
      {
        .young_selection_intensities(selection_proportion, rho_12)
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


  P_star <- .adjust_phenotypic_matrix(P, P1, b_M1, k1, tau, stage1_indices)
  C_star <- .adjust_genotypic_matrix(C, G1, b_M1, k1, tau, P1, stage1_indices)


  stage1_metrics <- .stage_metrics(b_M1, P1, G1, w1, k1)
  stage2_metrics <- .stage_metrics(b_M2, P_star, C_star, w, k2)


  gain_ratios_1 <- stage1_metrics$E / d1
  gain_ratios_1[!is.finite(gain_ratios_1)] <- NA_real_

  gain_ratios_2 <- stage2_metrics$E / d2
  gain_ratios_2[!is.finite(gain_ratios_2)] <- NA_real_


  trait_names_1 <- if (!is.null(colnames(P1))) colnames(P1) else paste0("Trait", 1:n1)
  trait_names_2 <- if (!is.null(colnames(P))) colnames(P) else paste0("Trait", 1:n)

  summary_stage1 <- data.frame(
    Stage = 1,
    Trait = trait_names_1,
    b_M = round(as.numeric(b_M1), 4),
    d = round(d1, 4),
    E = round(stage1_metrics$E, 4),
    Ratio = round(gain_ratios_1, 4),
    stringsAsFactors = FALSE
  )

  summary_stage2 <- data.frame(
    Stage = 2,
    Trait = trait_names_2,
    b_M = round(as.numeric(b_M2), 4),
    d = round(d2, 4),
    E = round(stage2_metrics$E, 4),
    Ratio = round(gain_ratios_2, 4),
    stringsAsFactors = FALSE
  )


  result <- list(
    b_M1 = as.numeric(b_M1),
    b_M2 = as.numeric(b_M2),
    b_R1 = as.numeric(b_R1),
    b_R2 = as.numeric(b_R2),
    b1 = as.numeric(b1),
    b2 = as.numeric(b2),
    K_M1 = K_M1,
    K_M2 = K_M2,
    theta1 = theta1,
    theta2 = theta2,
    d1 = d1,
    d2 = d2,
    gain_ratios_1 = setNames(gain_ratios_1, trait_names_1),
    gain_ratios_2 = setNames(gain_ratios_2, trait_names_2),
    stage1_metrics = list(
      R = stage1_metrics$R,
      E = setNames(stage1_metrics$E, trait_names_1),
      rho_H = stage1_metrics$rho_H,
      sigma_I = stage1_metrics$sigma_I
    ),
    stage2_metrics = list(
      R = stage2_metrics$R,
      E = setNames(stage2_metrics$E, trait_names_2),
      rho_H = stage2_metrics$rho_H,
      sigma_I = stage2_metrics$sigma_I
    ),
    P_star = P_star,
    C_star = C_star,
    rho_12 = rho_12,
    k1 = k1,
    k2 = k2,
    tau = tau,
    selection_proportion = selection_proportion,
    summary_stage1 = summary_stage1,
    summary_stage2 = summary_stage2
  )

  class(result) <- c("mppg_lpsi", "multistage_index", "list")

  result
}
