#' Genomic Selection Indices
#' @name genomic_indices
#'
#' @description
#' Implements genomic selection indices using GEBVs (Genomic Estimated Breeding Values).
#'
#' Methods included:
#' - Linear Genomic Selection Index (LGSI)
#' - Combined Linear Genomic Selection Index (CLGSI - phenotypes + GEBVs)
#'
#' @keywords internal
#' @importFrom stats cov
NULL

#' Linear Genomic Selection Index (LGSI)
#'
#' @description
#' Implements the Linear Genomic Selection Index where selection is based solely on
#' Genomic Estimated Breeding Values (GEBVs). This is used for selecting candidates
#' that have been genotyped but not phenotyped (e.g., in a testing population).
#'
#' @param gebv_mat Matrix of GEBVs (n_genotypes x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param wmat Economic weights matrix (n_traits x k), or vector
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param reliability Optional. Reliability of GEBVs (correlation between GEBV and true BV).
#'   Can be:
#'   - Single value (applied to all traits)
#'   - Vector of length n_traits (one per trait)
#'   - NULL (default): estimated from GEBV variance (assumes reliability = GEBV_var / G_var)
#' @param selection_intensity Selection intensity i (default: 2.063 for 10\% selection)
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{b} - Index coefficients
#'     \item \code{P_gebv} - GEBV variance-covariance matrix
#'     \item \code{reliability} - Reliability values used
#'     \item \code{Delta_H} - Expected genetic advance per trait
#'     \item \code{GA} - Overall genetic advance in the index
#'     \item \code{PRE} - Percent relative efficiency (if GAY provided)
#'     \item \code{hI2} - Index heritability
#'     \item \code{rHI} - Index accuracy
#'     \item \code{sigma_I} - Standard deviation of the index
#'     \item \code{summary} - Data frame with all metrics
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The LGSI maximizes the correlation between the index I = b' * gebv and
#' the aggregate genotype H = w' * g.
#'
#' Index coefficients: \eqn{\mathbf{b} = \mathbf{P}_{\hat{g}}^{-1} \mathbf{C}_{\hat{g}g} \mathbf{w}}
#'
#' Where:
#' - \eqn{\mathbf{P}_{\hat{g}}} = Var(gebv) - variance-covariance of GEBVs
#' - \eqn{\mathbf{C}_{\hat{g}g}} = Cov(gebv, g) - covariance between GEBVs and true breeding values
#'
#' If reliability (r) is known: \eqn{\mathbf{C}_{\hat{g}g} = \text{diag}(r) \mathbf{P}_{\hat{g}}}
#'
#' Expected response: \eqn{\Delta \mathbf{H} = \frac{i}{\sigma_I} \mathbf{C}_{\hat{g}g} \mathbf{b}}
#'
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing.
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
#'
#' # Simulate GEBVs (in practice, these come from genomic prediction)
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- ncol(gmat)
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
#'                    nrow = n_genotypes, ncol = n_traits)
#' colnames(gebv_mat) <- colnames(gmat)
#'
#' # Economic weights
#' weights <- c(10, 5, 3, 3, 5, 8, 4)
#'
#' # Calculate LGSI
#' result <- lgsi(gebv_mat, gmat, weights, reliability = 0.7)
#' print(result$summary)
#' }
lgsi <- function(gebv_mat, gmat, wmat, wcol = 1,
                 reliability = NULL,
                 selection_intensity = 2.063,
                 GAY = NULL) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  gebv_mat <- as.matrix(gebv_mat)
  gmat <- as.matrix(gmat)

  n_genotypes <- nrow(gebv_mat)
  n_traits <- ncol(gebv_mat)

  if (nrow(gmat) != n_traits || ncol(gmat) != n_traits) {
    stop("gmat dimensions must match number of traits in gebv_mat")
  }

  # Check for NAs
  if (any(is.na(gebv_mat))) {
    stop("gebv_mat contains NA values. Please remove or impute missing values.")
  }

  # Handle wmat
  if (is.vector(wmat)) {
    wmat <- matrix(wmat, ncol = 1)
  } else {
    wmat <- as.matrix(wmat)
  }

  if (nrow(wmat) != n_traits) {
    stop("Number of rows in wmat must equal number of traits")
  }

  if (wcol < 1 || wcol > ncol(wmat)) {
    stop("wcol must be between 1 and ", ncol(wmat))
  }

  w <- as.numeric(wmat[, wcol])

  # ============================================================================
  # STEP 1: Compute GEBV Variance-Covariance Matrix (using R's cov)
  # ============================================================================

  P_gebv <- cov(gebv_mat)

  # Check for zero variance traits
  diag_var <- diag(P_gebv)
  if (any(diag_var < 1e-10)) {
    zero_var_traits <- which(diag_var < 1e-10)
    warning("Trait(s) ", paste(zero_var_traits, collapse = ", "),
            " have near-zero GEBV variance")
  }

  # ============================================================================
  # STEP 2: Handle Reliability
  # ============================================================================

  if (is.null(reliability)) {
    # Estimate reliability as ratio of GEBV variance to genetic variance
    # This assumes: Var(GEBV) = r^2 * Var(G)
    # So: r = sqrt(Var(GEBV) / Var(G))
    diag_g <- diag(gmat)
    r_vec <- sqrt(pmin(diag_var / diag_g, 1))  # Clamp to [0, 1]

    # Warn if any reliabilities are very low
    if (any(r_vec < 0.3)) {
      warning("Estimated reliabilities are low for some traits. Consider providing known reliability values.")
    }
  } else if (length(reliability) == 1) {
    # Single reliability value for all traits
    if (reliability < 0 || reliability > 1) {
      stop("Reliability must be between 0 and 1")
    }
    r_vec <- rep(reliability, n_traits)
  } else if (length(reliability) == n_traits) {
    # Vector of reliabilities
    r_vec <- as.numeric(reliability)
    if (any(r_vec < 0 | r_vec > 1)) {
      stop("All reliability values must be between 0 and 1")
    }
  } else {
    stop("reliability must be NULL, a single value, or a vector of length n_traits")
  }

  # ============================================================================
  # STEP 3: Compute Covariance between GEBV and True Breeding Values
  # C_gebv_g = diag(r) * P_gebv (approximation based on reliability)
  # ============================================================================

  C_gebv_g <- P_gebv * outer(r_vec, r_vec, "*")
  # For diagonal: C_gebv_g = diag(r_vec) %*% P_gebv %*% diag(r_vec)
  # But simpler: element-wise multiplication with reliability matrix

  # Alternative: just scale by reliability
  # C_gebv_g <- sweep(sweep(P_gebv, 1, r_vec, "*"), 2, r_vec, "*")

  # ============================================================================
  # STEP 4: Solve for Index Coefficients
  # b = P_gebv^(-1) * C_gebv_g * w
  # ============================================================================

  # Use existing helper function
  P_inv_C <- .solve_sym_multi(P_gebv, C_gebv_g)
  b <- P_inv_C %*% w
  b <- as.numeric(b)

  # Add trait names to coefficients
  if (!is.null(colnames(gebv_mat))) {
    names(b) <- colnames(gebv_mat)
  }

  # Check for numerical issues
  if (any(is.na(b)) || any(is.infinite(b))) {
    stop("Index coefficients contain NA or Inf. Check that P_gebv is invertible.")
  }

  # ============================================================================
  # STEP 5: Calculate Index Metrics
  # ============================================================================

  # Use existing C++ primitives for calculations
  bPb <- cpp_quadratic_form_sym(b, P_gebv)
  bCb <- cpp_quadratic_form_sym(b, C_gebv_g)
  sigma_I <- if (bPb > 0) sqrt(bPb) else NA_real_

  # Index heritability
  hI2 <- if (!is.na(bPb) && bPb > 0) bCb / bPb else NA_real_

  # Index accuracy
  rHI <- if (!is.na(hI2) && hI2 >= 0 && hI2 <= 1) sqrt(hI2) else NA_real_

  # Expected genetic response per trait
  Delta_H_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    selection_intensity * (C_gebv_g %*% b) / sigma_I
  } else {
    rep(NA_real_, n_traits)
  }

  # Overall genetic advance
  bCw <- cpp_quadratic_form(b, C_gebv_g, w)
  GA <- if (!is.na(sigma_I) && sigma_I > 0) {
    selection_intensity * bCw / sigma_I
  } else {
    NA_real_
  }

  # Percent relative efficiency
  PRE <- if (!is.null(GAY) && !is.na(GA)) {
    (GA / GAY) * 100
  } else if (!is.na(GA)) {
    GA * 100
  } else {
    NA_real_
  }

  # ============================================================================
  # STEP 6: Build Summary Data Frame
  # ============================================================================

  b_vec <- round(b, 4)
  b_df <- as.data.frame(matrix(b_vec, nrow = 1))
  colnames(b_df) <- paste0("b.", seq_len(length(b_vec)))

  summary_df <- data.frame(
    b_df,
    GA = round(GA, 4),
    PRE = round(PRE, 4),
    hI2 = round(hI2, 4),
    rHI = round(rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Name Delta_H vector
  trait_names <- colnames(gebv_mat)
  if (!is.null(trait_names)) {
    names(Delta_H_vec) <- trait_names
    names(r_vec) <- trait_names
  }

  # ============================================================================
  # STEP 7: Return Results
  # ============================================================================

  result <- list(
    b = b_vec,
    P_gebv = P_gebv,
    C_gebv_g = C_gebv_g,
    reliability = r_vec,
    Delta_H = as.numeric(Delta_H_vec),
    GA = GA,
    PRE = PRE,
    hI2 = hI2,
    rHI = rHI,
    sigma_I = sigma_I,
    selection_intensity = selection_intensity,
    summary = summary_df
  )

  class(result) <- c("lgsi", "selection_index", "list")

  return(result)
}


#' Combined Linear Genomic Selection Index (CLGSI)
#'
#' @description
#' Implements the Combined Linear Genomic Selection Index where selection combines
#' both phenotypic observations and Genomic Estimated Breeding Values (GEBVs).
#' This is used for selecting candidates with both phenotype and genotype data
#' (e.g., in a training population).
#'
#' @param phen_mat Matrix of phenotypes (n_genotypes x n_traits)
#' @param gebv_mat Matrix of GEBVs (n_genotypes x n_traits)
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits)
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param wmat Economic weights matrix (n_traits x k), or vector
#' @param wcol Weight column to use if wmat has multiple columns (default: 1)
#' @param reliability Optional. Reliability of GEBVs (see lgsi() for details)
#' @param selection_intensity Selection intensity i (default: 2.063 for 10\% selection)
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation
#'
#' @return List with components:
#'   \itemize{
#'     \item \code{b_y} - Coefficients for phenotypes
#'     \item \code{b_g} - Coefficients for GEBVs
#'     \item \code{b_combined} - Full coefficient vector [b_y; b_g]
#'     \item \code{P_combined} - Combined variance matrix
#'     \item \code{Delta_H} - Expected genetic advance per trait
#'     \item \code{GA} - Overall genetic advance
#'     \item \code{PRE} - Percent relative efficiency
#'     \item \code{hI2} - Index heritability
#'     \item \code{rHI} - Index accuracy
#'     \item \code{summary} - Data frame with all metrics
#'   }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The CLGSI combines phenotypic and genomic information:
#' \deqn{I = \mathbf{b}_y' \mathbf{y} + \mathbf{b}_g' \hat{\mathbf{g}}}
#'
#' Coefficients are obtained by solving the partitioned system:
#' \deqn{\begin{bmatrix} \mathbf{b}_y \\ \mathbf{b}_g \end{bmatrix} =
#'       \begin{bmatrix} \mathbf{P} & \mathbf{P}_{y\hat{g}} \\
#'                       \mathbf{P}_{y\hat{g}}' & \mathbf{P}_{\hat{g}} \end{bmatrix}^{-1}
#'       \begin{bmatrix} \mathbf{G} \\ \mathbf{C}_{\hat{g}g} \end{bmatrix} \mathbf{w}}
#'
#' Where:
#' - \eqn{\mathbf{P}} = Var(phenotypes)
#' - \eqn{\mathbf{P}_{\hat{g}}} = Var(GEBVs)
#' - \eqn{\mathbf{P}_{y\hat{g}}} = Cov(phenotypes, GEBVs)
#' - \eqn{\mathbf{G}} = Genotypic variance-covariance
#' - \eqn{\mathbf{C}_{\hat{g}g}} = Cov(GEBV, true BV)
#'
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern Plant Breeding.
#' Springer International Publishing.
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
#' n_traits <- ncol(gmat)
#'
#' phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
#'                    nrow = n_genotypes, ncol = n_traits)
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
#'                    nrow = n_genotypes, ncol = n_traits)
#' colnames(phen_mat) <- colnames(gmat)
#' colnames(gebv_mat) <- colnames(gmat)
#'
#' # Economic weights
#' weights <- c(10, 5, 3, 3, 5, 8, 4)
#'
#' # Calculate CLGSI
#' result <- clgsi(phen_mat, gebv_mat, pmat, gmat, weights, reliability = 0.7)
#' print(result$summary)
#' }
clgsi <- function(phen_mat, gebv_mat, pmat, gmat, wmat, wcol = 1,
                  reliability = NULL,
                  selection_intensity = 2.063,
                  GAY = NULL) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  phen_mat <- as.matrix(phen_mat)
  gebv_mat <- as.matrix(gebv_mat)
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)

  n_genotypes <- nrow(phen_mat)
  n_traits <- ncol(phen_mat)

  if (nrow(gebv_mat) != n_genotypes || ncol(gebv_mat) != n_traits) {
    stop("phen_mat and gebv_mat must have same dimensions")
  }

  if (nrow(pmat) != n_traits || ncol(pmat) != n_traits) {
    stop("pmat dimensions must match number of traits")
  }

  if (nrow(gmat) != n_traits || ncol(gmat) != n_traits) {
    stop("gmat dimensions must match number of traits")
  }

  # Check for NAs
  if (any(is.na(phen_mat))) {
    stop("phen_mat contains NA values")
  }
  if (any(is.na(gebv_mat))) {
    stop("gebv_mat contains NA values")
  }

  # Handle wmat
  if (is.vector(wmat)) {
    wmat <- matrix(wmat, ncol = 1)
  } else {
    wmat <- as.matrix(wmat)
  }

  if (nrow(wmat) != n_traits) {
    stop("Number of rows in wmat must equal number of traits")
  }

  if (wcol < 1 || wcol > ncol(wmat)) {
    stop("wcol must be between 1 and ", ncol(wmat))
  }

  w <- as.numeric(wmat[, wcol])

  # ============================================================================
  # STEP 1: Compute Variance-Covariance Matrices
  # ============================================================================

  # Use R's cov for all variance-covariance calculations
  P_gebv <- cov(gebv_mat)
  P_yg <- cov(phen_mat, gebv_mat)  # Covariance between phenotypes and GEBVs

  # ============================================================================
  # STEP 2: Handle Reliability
  # ============================================================================

  if (is.null(reliability)) {
    diag_var_gebv <- diag(P_gebv)
    diag_var_g <- diag(gmat)
    r_vec <- sqrt(pmin(diag_var_gebv / diag_var_g, 1))

    if (any(r_vec < 0.3)) {
      warning("Estimated reliabilities are low for some traits")
    }
  } else if (length(reliability) == 1) {
    if (reliability < 0 || reliability > 1) {
      stop("Reliability must be between 0 and 1")
    }
    r_vec <- rep(reliability, n_traits)
  } else if (length(reliability) == n_traits) {
    r_vec <- as.numeric(reliability)
    if (any(r_vec < 0 | r_vec > 1)) {
      stop("All reliability values must be between 0 and 1")
    }
  } else {
    stop("reliability must be NULL, a single value, or a vector of length n_traits")
  }

  # Covariance between GEBV and true BV
  C_gebv_g <- P_gebv * outer(r_vec, r_vec, "*")

  # ============================================================================
  # STEP 3: Build Combined 2t x 2t Variance Matrix
  # P_combined = [P,     P_yg  ]
  #              [P_yg', P_gebv]
  # ============================================================================

  P_combined <- rbind(
    cbind(pmat, P_yg),
    cbind(t(P_yg), P_gebv)
  )

  # Check symmetry
  if (max(abs(P_combined - t(P_combined))) > 1e-6) {
    warning("Combined variance matrix is not symmetric")
  }

  # ============================================================================
  # STEP 4: Build Combined 2t x t Genetic Covariance Matrix
  # G_combined = [G      ]
  #              [C_gebv_g]
  # ============================================================================

  G_combined <- rbind(gmat, C_gebv_g)

  # ============================================================================
  # STEP 5: Solve for Combined Index Coefficients
  # [b_y; b_g] = P_combined^(-1) * G_combined * w
  # ============================================================================

  # Compute RHS
  rhs <- G_combined %*% w

  # Solve using existing solve (which uses cpp_symmetric_solve internally)
  # For 2t x 2t system, we can iterate over columns of rhs
  b_combined <- .solve_sym_multi(P_combined, matrix(rhs, ncol = 1))
  b_combined <- as.numeric(b_combined)

  # Check for numerical issues
  if (any(is.na(b_combined)) || any(is.infinite(b_combined))) {
    stop("Index coefficients contain NA or Inf. Check matrix conditioning.")
  }

  # Split into phenotype and GEBV coefficients
  b_y <- b_combined[1:n_traits]
  b_g <- b_combined[(n_traits + 1):(2 * n_traits)]

  # ============================================================================
  # STEP 6: Calculate Index Metrics
  # ============================================================================

  # Index variance: b' * P_combined * b
  bPb <- cpp_quadratic_form_sym(b_combined, P_combined)
  sigma_I <- if (bPb > 0) sqrt(bPb) else NA_real_

  # Index-aggregate genotype covariance: b' * G_combined * w
  bGw <- as.numeric(t(b_combined) %*% G_combined %*% w)

  # Index heritability (approximation using combined system)
  bGb <- as.numeric(t(b_combined) %*% G_combined %*% b_combined[1:n_traits] +
                    t(b_combined) %*% rbind(matrix(0, n_traits, n_traits), C_gebv_g) %*% b_combined[(n_traits+1):(2*n_traits)])
  hI2 <- if (!is.na(bPb) && bPb > 0) bGb / bPb else NA_real_
  hI2 <- pmin(pmax(hI2, 0), 1)  # Clamp to [0, 1]

  # Index accuracy
  rHI <- if (!is.na(hI2) && hI2 >= 0) sqrt(hI2) else NA_real_

  # Expected genetic response per trait
  # Delta_H = (selection_intensity / sigma_I) * (G' * b_y + C_gebv_g' * b_g)
  Delta_H_vec <- if (!is.na(sigma_I) && sigma_I > 0) {
    selection_intensity * as.vector(t(G_combined) %*% b_combined) / sigma_I
  } else {
    rep(NA_real_, n_traits)
  }

  # Overall genetic advance
  GA <- if (!is.na(sigma_I) && sigma_I > 0) {
    selection_intensity * bGw / sigma_I
  } else {
    NA_real_
  }

  # Percent relative efficiency
  PRE <- if (!is.null(GAY) && !is.na(GA)) {
    (GA / GAY) * 100
  } else if (!is.na(GA)) {
    GA * 100
  } else {
    NA_real_
  }

  # ============================================================================
  # STEP 7: Build Summary Data Frame
  # ============================================================================

  # Coefficients for phenotypes
  b_y_vec <- round(b_y, 4)
  b_y_df <- as.data.frame(matrix(b_y_vec, nrow = 1))
  colnames(b_y_df) <- paste0("b_y.", seq_len(length(b_y_vec)))

  # Coefficients for GEBVs
  b_g_vec <- round(b_g, 4)
  b_g_df <- as.data.frame(matrix(b_g_vec, nrow = 1))
  colnames(b_g_df) <- paste0("b_g.", seq_len(length(b_g_vec)))

  summary_df <- data.frame(
    b_y_df,
    b_g_df,
    GA = round(GA, 4),
    PRE = round(PRE, 4),
    hI2 = round(hI2, 4),
    rHI = round(rHI, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Name vectors
  trait_names <- colnames(phen_mat)
  if (!is.null(trait_names)) {
    names(Delta_H_vec) <- trait_names
    names(r_vec) <- trait_names
    names(b_y) <- trait_names
    names(b_g) <- trait_names
  }

  # ============================================================================
  # STEP 8: Return Results
  # ============================================================================

  result <- list(
    b_y = b_y,
    b_g = b_g,
    b_combined = b_combined,
    P_combined = P_combined,
    G_combined = G_combined,
    reliability = r_vec,
    Delta_H = as.numeric(Delta_H_vec),
    GA = GA,
    PRE = PRE,
    hI2 = hI2,
    rHI = rHI,
    sigma_I = sigma_I,
    selection_intensity = selection_intensity,
    summary = summary_df
  )

  class(result) <- c("clgsi", "selection_index", "list")

  return(result)
}
