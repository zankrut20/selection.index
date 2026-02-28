#' Linear Marker and Genome-Wide Selection Indices (Chapter 4)
#' @name marker_indices
#'
#' @description
#' Implements marker-based selection indices from Chapter 4 (Lande & Thompson, 1990).
#' These methods incorporate molecular marker information directly into selection indices.
#'
#' Methods included:
#' - LMSI    : Linear Marker Selection Index (Section 4.1)
#' - GW-LMSI : Genome-Wide Linear Marker Selection Index (Section 4.2)
#'
#' @keywords internal
#' @importFrom stats cov
#' @importFrom MASS ginv
NULL


#' Linear Marker Selection Index (LMSI)
#'
#' @description
#' Implements the LMSI which combines phenotypic information with molecular
#' marker scores from statistically significant markers (Lande & Thompson, 1990).
#' The index is I = b_y' * y + b_s' * s, where y are phenotypes and s are
#' marker scores.
#'
#' @param phen_mat Matrix of phenotypes (n_genotypes x n_traits).
#'   Can be NULL if G_s is provided directly (theoretical case where
#'   covariance structure is known without needing empirical data).
#' @param marker_scores Matrix of marker scores (n_genotypes x n_traits).
#'   These are computed as s_j = sum(x_jk * beta_jk) where x_jk is the coded
#'   marker value and beta_jk is the estimated marker effect for trait j.
#'   Can be NULL if G_s is provided directly.
#' @param pmat Phenotypic variance-covariance matrix (n_traits x n_traits).
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param G_s Genomic covariance matrix explained by markers (n_traits x n_traits).
#'   This represents Var(s) which approximates Cov(y, s) when markers fully
#'   explain genetic variance. If provided, phen_mat and marker_scores become
#'   optional as the covariance structure is specified directly.
#'   If NULL, computed empirically from marker_scores and phen_mat.
#' @param wmat Economic weights matrix (n_traits x k), or vector.
#' @param wcol Weight column to use if wmat has multiple columns (default: 1).
#' @param selection_intensity Selection intensity k (default: 2.063 for 10\% selection).
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation.
#'
#' @return List of class \code{"lmsi"} with components:
#' \describe{
#'   \item{\code{b_y}}{Coefficients for phenotypes (n_traits vector).}
#'   \item{\code{b_s}}{Coefficients for marker scores (n_traits vector).}
#'   \item{\code{b_combined}}{Combined coefficient vector [b_y; b_s] (2*n_traits vector).}
#'   \item{\code{P_L}}{Combined phenotypic-marker covariance matrix (2*n_traits x 2*n_traits).}
#'   \item{\code{G_L}}{Combined genetic-marker covariance matrix (2*n_traits x n_traits).}
#'   \item{\code{G_s}}{Genomic covariance matrix explained by markers (n_traits x n_traits).}
#'   \item{\code{rHI}}{Index accuracy (correlation between index and breeding objective).}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{R}}{Selection response (k * sigma_I).}
#'   \item{\code{Delta_H}}{Expected genetic gain per trait (vector of length n_traits).}
#'   \item{\code{GA}}{Overall genetic advance in breeding objective.}
#'   \item{\code{PRE}}{Percent relative efficiency (if GAY provided).}
#'   \item{\code{hI2}}{Index heritability.}
#'   \item{\code{summary}}{Data frame with coefficients and metrics (combined view).}
#'   \item{\code{phenotype_coeffs}}{Data frame with phenotype coefficients only.}
#'   \item{\code{marker_coeffs}}{Data frame with marker score coefficients only.}
#'   \item{\code{coeff_analysis}}{Data frame with coefficient distribution analysis.}
#' }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The LMSI maximizes the correlation between the index
#' \eqn{I_{LMSI} = \mathbf{b}_y^{\prime}\mathbf{y} + \mathbf{b}_s^{\prime}\mathbf{s}}
#' and the breeding objective \eqn{H = \mathbf{w}^{\prime}\mathbf{g}}.
#'
#' Combined covariance matrices:
#' \deqn{\mathbf{P}_L = \begin{bmatrix} \mathbf{P} & \text{Cov}(\mathbf{y}, \mathbf{s}) \\ \text{Cov}(\mathbf{y}, \mathbf{s})^{\prime} & \text{Var}(\mathbf{s}) \end{bmatrix}}
#' \deqn{\mathbf{G}_L = \begin{bmatrix} \mathbf{G} \\ \mathbf{G}_s \end{bmatrix}}
#'
#' where \eqn{\mathbf{P}} is the phenotypic variance, \eqn{\text{Cov}(\mathbf{y}, \mathbf{s})}
#' is the covariance between phenotypes and marker scores (computed from data),
#' \eqn{\text{Var}(\mathbf{s})} is the variance of marker scores, \eqn{\mathbf{G}} is the
#' genotypic variance, and \eqn{\mathbf{G}_s} represents the genetic covariance
#' explained by markers.
#'
#' Index coefficients:
#' \deqn{\mathbf{b}_{LMSI} = \mathbf{P}_L^{-1} \mathbf{G}_L \mathbf{w}}
#'
#' Accuracy:
#' \deqn{\rho_{HI} = \sqrt{\frac{\mathbf{b}_{LMSI}^{\prime} \mathbf{G}_L \mathbf{w}}{\mathbf{w}^{\prime} \mathbf{G} \mathbf{w}}}}
#'
#' Selection response:
#' \deqn{R_{LMSI} = k \sigma_{I_{LMSI}} = k \sqrt{\mathbf{b}_{LMSI}^{\prime} \mathbf{P}_L \mathbf{b}_{LMSI}}}
#'
#' Expected genetic gain per trait:
#' \deqn{\mathbf{E}_{LMSI} = k \frac{\mathbf{G}_L^{\prime} \mathbf{b}_{LMSI}}{\sigma_{I_{LMSI}}}}
#'
#' @references
#' Lande, R., & Thompson, R. (1990). Efficiency of marker-assisted selection
#' in the improvement of quantitative traits. Genetics, 124(3), 743-756.
#'
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Chapter 4.
#'
#' @export
#' @examples
#' \dontrun{
#' # Load data
#' data(seldata)
#' pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate marker scores (in practice, computed from QTL mapping)
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- ncol(gmat)
#' marker_scores <- matrix(rnorm(n_genotypes * n_traits, mean = 5, sd = 1.5),
#'   nrow = n_genotypes, ncol = n_traits
#' )
#' colnames(marker_scores) <- colnames(gmat)
#'
#' # Simulate phenotypes
#' phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
#'   nrow = n_genotypes, ncol = n_traits
#' )
#' colnames(phen_mat) <- colnames(gmat)
#'
#' # Economic weights
#' weights <- c(10, 5, 3, 3, 5, 8, 4)
#'
#' # Calculate LMSI
#' result <- lmsi(phen_mat, marker_scores, pmat, gmat,
#'   G_s = NULL, wmat = weights
#' )
#' print(result$summary)
#' }
lmsi <- function(phen_mat = NULL, marker_scores = NULL,
                 pmat, gmat, G_s = NULL,
                 wmat, wcol = 1,
                 selection_intensity = 2.063,
                 GAY = NULL) {
  pmat <- as.matrix(pmat)
  gmat <- as.matrix(gmat)

  n_traits <- nrow(pmat)

  if (nrow(pmat) != ncol(pmat) || nrow(gmat) != ncol(gmat)) {
    stop("pmat and gmat must be square matrices")
  }

  if (nrow(pmat) != nrow(gmat)) {
    stop("pmat and gmat must have the same dimensions")
  }

  if (is.matrix(wmat)) {
    if (wcol > ncol(wmat)) {
      stop("wcol exceeds number of columns in wmat")
    }
    w <- as.vector(wmat[, wcol])
  } else {
    w <- as.vector(wmat)
  }

  if (length(w) != n_traits) {
    stop("Length of weights must equal number of traits")
  }

  if (is.null(G_s)) {
    if (is.null(phen_mat) || is.null(marker_scores)) {
      stop("Either G_s must be provided, or both phen_mat and marker_scores must be provided to compute covariance matrices")
    }

    phen_mat <- as.matrix(phen_mat)
    marker_scores <- as.matrix(marker_scores)

    if (ncol(phen_mat) != n_traits || ncol(marker_scores) != n_traits) {
      stop("Number of columns in phen_mat and marker_scores must equal n_traits")
    }

    if (nrow(phen_mat) != nrow(marker_scores)) {
      stop("phen_mat and marker_scores must have the same number of rows")
    }

    Cov_ys <- cov(phen_mat, marker_scores)

    Var_s <- cov(marker_scores)

    G_s <- Cov_ys
  } else {
    G_s <- as.matrix(G_s)
    if (nrow(G_s) != n_traits || ncol(G_s) != n_traits) {
      stop("G_s must be n_traits x n_traits matrix")
    }

    Cov_ys <- G_s
    Var_s <- G_s
  }


  P_L <- rbind(
    cbind(pmat, Cov_ys),
    cbind(t(Cov_ys), Var_s)
  )

  G_L <- rbind(gmat, G_s)


  G_L_w <- G_L %*% w

  b_combined <- tryCatch(
    {
      solve(P_L, G_L_w)
    },
    error = function(e) {
      warning("P_L is singular or near-singular, using generalized inverse")
      MASS::ginv(P_L) %*% G_L_w
    }
  )

  b_combined <- as.vector(b_combined)

  b_y <- b_combined[1:n_traits]
  b_s <- b_combined[(n_traits + 1):(2 * n_traits)]


  sigma_I_sq <- cpp_quadratic_form_sym(b_combined, P_L)
  sigma_I <- sqrt(max(sigma_I_sq, 0))

  numerator <- cpp_quadratic_form(b_combined, G_L, w)

  denominator <- cpp_quadratic_form_sym(w, gmat)


  rHI <- if (denominator > 0) {
    ratio <- max(0, min(numerator / denominator, 1.0)) # Cap at [0, 1]
    sqrt(ratio)
  } else {
    0
  }

  R <- selection_intensity * sigma_I

  if (sigma_I > 0) {
    Delta_H <- (selection_intensity / sigma_I) * as.vector(t(G_L) %*% b_combined)
  } else {
    Delta_H <- rep(0, n_traits)
  }

  GA <- sum(w * Delta_H)

  PRE <- if (!is.null(GAY) && !is.na(GAY) && GAY != 0) {
    (GA / GAY) * 100
  } else {
    NA_real_
  }

  hI2 <- if (sigma_I_sq > 0) min(numerator / sigma_I_sq, 1.0) else 0


  trait_names <- colnames(pmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait", 1:n_traits)
  }

  summary_df <- data.frame(
    Trait = rep(trait_names, 2),
    Component = rep(c("Phenotype", "MarkerScore"), each = n_traits),
    b = c(b_y, b_s),
    w = rep(w, 2),
    Delta_H = rep(Delta_H, 2),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  phenotype_summary <- data.frame(
    Trait = trait_names,
    b_phenotype = round(b_y, 6),
    weight = w,
    Delta_H = round(Delta_H, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  marker_summary <- data.frame(
    Trait = trait_names,
    b_marker = round(b_s, 6),
    weight = w,
    Delta_H = round(Delta_H, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  coeff_analysis <- data.frame(
    Component = c("Phenotype", "MarkerScore", "Combined"),
    Sum_Abs_Coeff = c(sum(abs(b_y)), sum(abs(b_s)), sum(abs(b_combined))),
    Mean_Abs_Coeff = c(mean(abs(b_y)), mean(abs(b_s)), mean(abs(b_combined))),
    Max_Abs_Coeff = c(max(abs(b_y)), max(abs(b_s)), max(abs(b_combined))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  attr(summary_df, "metrics") <- data.frame(
    rHI = round(rHI, 4),
    hI2 = round(hI2, 4),
    sigma_I = round(sigma_I, 4),
    R = round(R, 4),
    GA = round(GA, 4),
    PRE = if (!is.na(PRE)) round(PRE, 2) else NA_real_,
    stringsAsFactors = FALSE
  )


  result <- list(
    b_y = b_y, # Phenotype coefficients
    b_s = b_s, # Marker score coefficients
    b_combined = b_combined, # Full combined vector [b_y; b_s]

    P_L = P_L, # Combined covariance matrix
    G_L = G_L, # Combined genetic covariance matrix
    G_s = G_s, # Genetic covariance explained by markers

    rHI = rHI, # Accuracy
    hI2 = hI2, # Heritability
    sigma_I = sigma_I, # Index standard deviation
    R = R, # Selection response
    GA = GA, # Overall genetic advance
    PRE = PRE, # Percent relative efficiency

    Delta_H = Delta_H, # Expected genetic gain per trait

    selection_intensity = selection_intensity,
    trait_names = trait_names,
    summary = summary_df, # Combined summary
    phenotype_coeffs = phenotype_summary, # Phenotype coefficients only
    marker_coeffs = marker_summary, # Marker coefficients only
    coeff_analysis = coeff_analysis # Coefficient distribution analysis
  )

  class(result) <- c("lmsi", "marker_index", "list")
  result
}


#' Genome-Wide Linear Marker Selection Index (GW-LMSI)
#'
#' @description
#' Implements the GW-LMSI which uses all available genome-wide markers directly
#' as predictors in the selection index. Unlike LMSI which uses aggregated marker
#' scores per trait, GW-LMSI treats each individual marker as a separate predictor.
#'
#' @param marker_mat Matrix of marker genotypes (n_genotypes x n_markers).
#'   Typically coded as -1, 0, 1 or 0, 1, 2.
#' @param trait_mat Matrix of trait values (n_genotypes x n_traits).
#'   Used to compute G_GW if not provided.
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits).
#' @param P_GW Marker covariance matrix (n_markers x n_markers).
#'   If NULL, computed as Var(marker_mat).
#' @param G_GW Covariance between markers and traits (n_markers x n_traits).
#'   If NULL, computed as Cov(marker_mat, trait_mat).
#' @param wmat Economic weights matrix (n_traits x k), or vector.
#' @param wcol Weight column to use if wmat has multiple columns (default: 1).
#' @param lambda Ridge regularization parameter (default: 0).
#'   If lambda > 0, uses P_GW + lambda*I for regularization.
#'   Automatic warnings issued when n_markers > n_genotypes (high-dimensional case)
#'   or when P_GW is ill-conditioned. Recommended values: 0.01-0.1 times mean(diag(P_GW)).
#' @param selection_intensity Selection intensity k (default: 2.063 for 10\% selection).
#' @param GAY Optional. Genetic advance of comparative trait for PRE calculation.
#'
#' @return List of class \code{"gw_lmsi"} with components:
#' \describe{
#'   \item{\code{b}}{Index coefficients for markers (n_markers).}
#'   \item{\code{P_GW}}{Marker covariance matrix (n_markers x n_markers).}
#'   \item{\code{G_GW}}{Covariance between markers and traits (n_markers x n_traits).}
#'   \item{\code{rHI}}{Index accuracy (correlation between index and breeding objective).}
#'   \item{\code{sigma_I}}{Standard deviation of the index.}
#'   \item{\code{R}}{Selection response (k * sigma_I).}
#'   \item{\code{Delta_H}}{Expected genetic gain per trait (vector of length n_traits).}
#'   \item{\code{GA}}{Overall genetic advance in breeding objective.}
#'   \item{\code{PRE}}{Percent relative efficiency (if GAY provided).}
#'   \item{\code{hI2}}{Index heritability.}
#'   \item{\code{lambda}}{Ridge regularization parameter used.}
#'   \item{\code{n_markers}}{Number of markers.}
#'   \item{\code{high_dimensional}}{Logical indicating if n_markers > n_genotypes.}
#'   \item{\code{condition_number}}{Condition number of P_GW (if computed).}
#'   \item{\code{summary}}{Data frame with metrics.}
#' }
#'
#' @details
#' \strong{Mathematical Formulation:}
#'
#' The GW-LMSI maximizes the correlation between the index
#' \eqn{I_{GW} = \mathbf{b}_{GW}^{\prime}\mathbf{m}} and the breeding
#' objective \eqn{H = \mathbf{w}^{\prime}\mathbf{g}}.
#'
#' Marker covariance matrix:
#' \deqn{\mathbf{P}_{GW} = Var(\mathbf{m})}
#'
#' Covariance between markers and traits:
#' \deqn{\mathbf{G}_{GW} = Cov(\mathbf{m}, \mathbf{g})}
#'
#' Index coefficients (with optional Ridge regularization):
#' \deqn{\mathbf{b}_{GW} = (\mathbf{P}_{GW} + \lambda \mathbf{I})^{-1} \mathbf{G}_{GW} \mathbf{w}}
#'
#' Accuracy:
#' \deqn{\rho_{HI} = \sqrt{\frac{\mathbf{b}_{GW}^{\prime} \mathbf{G}_{GW} \mathbf{w}}{\mathbf{w}^{\prime} \mathbf{G} \mathbf{w}}}}
#'
#' Selection response:
#' \deqn{R_{GW} = k \sqrt{\mathbf{b}_{GW}^{\prime} \mathbf{P}_{GW} \mathbf{b}_{GW}}}
#'
#' Expected genetic gain per trait:
#' \deqn{\mathbf{E}_{GW} = k \frac{\mathbf{G}_{GW}^{\prime} \mathbf{b}_{GW}}{\sigma_{I_{GW}}}}
#'
#' \strong{Note on Singularity Detection and Regularization:}
#' The function automatically detects problematic cases:
#'
#' 1. **High-dimensional case**: When n_markers > n_genotypes, P_GW is mathematically
#'    singular (rank-deficient). The function issues a warning and suggests an
#'    appropriate lambda value.
#'
#' 2. **Ill-conditioned case**: When P_GW has a high condition number (> 1e10),
#'    indicating numerical instability.
#'
#' 3. **Numerical singularity**: When P_GW has eigenvalues near zero.
#'
#' Ridge regularization adds \eqn{\lambda}I to P_GW, ensuring positive definiteness. Recommended
#' lambda values are 0.01-0.1 times the average diagonal element of P_GW. Users can
#' also set lambda = 0 to force generalized inverse (less stable but sometimes needed).
#'
#' @references
#' Lande, R., & Thompson, R. (1990). Efficiency of marker-assisted selection
#' in the improvement of quantitative traits. Genetics, 124(3), 743-756.
#'
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Chapter 4.
#'
#' @export
#' @examples
#' \dontrun{
#' # Load data
#' data(seldata)
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate marker data
#' set.seed(123)
#' n_genotypes <- 100
#' n_markers <- 200
#' n_traits <- ncol(gmat)
#'
#' # Marker matrix (coded as 0, 1, 2)
#' marker_mat <- matrix(sample(0:2, n_genotypes * n_markers, replace = TRUE),
#'   nrow = n_genotypes, ncol = n_markers
#' )
#'
#' # Trait matrix
#' trait_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
#'   nrow = n_genotypes, ncol = n_traits
#' )
#'
#' # Economic weights
#' weights <- c(10, 5, 3, 3, 5, 8, 4)
#'
#' # Calculate GW-LMSI with Ridge regularization
#' result <- gw_lmsi(marker_mat, trait_mat, gmat,
#'   wmat = weights, lambda = 0.01
#' )
#' print(result$summary)
#' }
gw_lmsi <- function(marker_mat, trait_mat = NULL,
                    gmat,
                    P_GW = NULL, G_GW = NULL,
                    wmat, wcol = 1,
                    lambda = 0,
                    selection_intensity = 2.063,
                    GAY = NULL) {
  marker_mat <- as.matrix(marker_mat)
  gmat <- as.matrix(gmat)

  n_genotypes <- nrow(marker_mat)
  n_markers <- ncol(marker_mat)
  n_traits <- nrow(gmat)

  if (nrow(gmat) != ncol(gmat)) {
    stop("gmat must be a square matrix")
  }

  if (is.matrix(wmat)) {
    if (wcol > ncol(wmat)) {
      stop("wcol exceeds number of columns in wmat")
    }
    w <- as.vector(wmat[, wcol])
  } else {
    w <- as.vector(wmat)
  }

  if (length(w) != n_traits) {
    stop("Length of weights must equal number of traits")
  }

  if (is.null(P_GW) || is.null(G_GW)) {
    if (is.null(trait_mat)) {
      stop("Either (P_GW and G_GW) or trait_mat must be provided")
    }

    trait_mat <- as.matrix(trait_mat)

    if (ncol(trait_mat) != n_traits) {
      stop("Number of columns in trait_mat must equal n_traits")
    }

    if (nrow(trait_mat) != n_genotypes) {
      stop("marker_mat and trait_mat must have the same number of rows")
    }

    if (is.null(P_GW)) {
      P_GW <- cov(marker_mat)
    }

    if (is.null(G_GW)) {
      G_GW <- cov(marker_mat, trait_mat)
    }
  } else {
    P_GW <- as.matrix(P_GW)
    G_GW <- as.matrix(G_GW)

    if (nrow(P_GW) != n_markers || ncol(P_GW) != n_markers) {
      stop("P_GW must be n_markers x n_markers matrix")
    }

    if (nrow(G_GW) != n_markers || ncol(G_GW) != n_traits) {
      stop("G_GW must be n_markers x n_traits matrix")
    }
  }

  if (lambda < 0) {
    stop("lambda must be non-negative")
  }


  high_dimensional <- (n_markers > n_genotypes)

  condition_number <- NA_real_

  tryCatch(
    {
      P_GW_eigs <- eigen(P_GW, symmetric = TRUE, only.values = TRUE)$values
      max_eig <- max(P_GW_eigs)
      min_eig_pos <- P_GW_eigs[P_GW_eigs > 1e-14] # Exclude numerical zeros

      if (length(min_eig_pos) > 0) {
        min_eig <- min(min_eig_pos)
        condition_number <- max_eig / min_eig
      }
    },
    error = function(e) {
      condition_number <<- NA_real_
    }
  )

  if (high_dimensional && lambda == 0) {
    avg_diag <- mean(diag(P_GW))
    suggested_lambda <- 0.01 * avg_diag

    warning(
      sprintf(
        "High-dimensional case detected (n_markers = %d > n_genotypes = %d).\n  P_GW is singular. Ridge regularization is required.\n  Consider setting lambda > 0 (suggested: lambda = %.4f).\n  Currently using generalized inverse, which may be unstable.",
        n_markers, n_genotypes, suggested_lambda
      )
    )
  }

  if (!is.na(condition_number)) {
    if (condition_number > 1e10) {
      warning(
        sprintf(
          "P_GW is ill-conditioned (condition number = %.2e).\n  Consider using lambda > 0 for numerical stability.",
          condition_number
        )
      )
    }
  } else if (!high_dimensional && lambda == 0) {
    warning("P_GW appears to be numerically singular. Consider using lambda > 0 for stability.")
  }


  if (lambda > 0) {
    P_GW_reg <- P_GW + lambda * diag(n_markers)
    ridge_applied <- TRUE
  } else {
    P_GW_reg <- P_GW
    ridge_applied <- FALSE
  }


  G_GW_w <- G_GW %*% w

  b <- tryCatch(
    {
      solve(P_GW_reg, G_GW_w)
    },
    error = function(e) {
      if (!high_dimensional) {
        warning("P_GW is singular or near-singular, using generalized inverse")
      }
      if (!ridge_applied && !high_dimensional) {
        warning("Consider using lambda > 0 for Ridge regularization")
      }
      MASS::ginv(P_GW_reg) %*% G_GW_w
    }
  )

  b <- as.vector(b)


  sigma_I_sq <- cpp_quadratic_form_sym(b, P_GW)
  sigma_I <- sqrt(max(sigma_I_sq, 0))

  numerator <- cpp_quadratic_form(b, G_GW, w)

  denominator <- cpp_quadratic_form_sym(w, gmat)


  rHI <- if (denominator > 0) {
    ratio <- max(0, min(numerator / denominator, 1.0)) # Cap at [0, 1]
    sqrt(ratio)
  } else {
    0
  }

  R <- selection_intensity * sigma_I

  if (sigma_I > 0) {
    Delta_H <- (selection_intensity / sigma_I) * as.vector(t(G_GW) %*% b)
  } else {
    Delta_H <- rep(0, n_traits)
  }

  GA <- sum(w * Delta_H)

  PRE <- if (!is.null(GAY) && !is.na(GAY) && GAY != 0) {
    (GA / GAY) * 100
  } else {
    NA_real_
  }

  hI2 <- if (sigma_I_sq > 0) min(numerator / sigma_I_sq, 1.0) else 0


  trait_names <- colnames(gmat)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait", 1:n_traits)
  }

  summary_df <- data.frame(
    Metric = c("rHI", "hI2", "sigma_I", "R", "GA", "PRE", "n_markers", "lambda"),
    Value = c(
      round(rHI, 4),
      round(hI2, 4),
      round(sigma_I, 4),
      round(R, 4),
      round(GA, 4),
      if (!is.na(PRE)) round(PRE, 2) else NA_real_,
      n_markers,
      lambda
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  trait_gains <- data.frame(
    Trait = trait_names,
    w = w,
    Delta_H = round(Delta_H, 4),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )


  result <- list(
    b = b,
    P_GW = P_GW,
    G_GW = G_GW,
    rHI = rHI,
    hI2 = hI2,
    sigma_I = sigma_I,
    R = R,
    Delta_H = Delta_H,
    GA = GA,
    PRE = PRE,
    lambda = lambda,
    ridge_applied = ridge_applied,
    high_dimensional = high_dimensional,
    condition_number = condition_number,
    n_markers = n_markers,
    n_genotypes = n_genotypes,
    n_traits = n_traits,
    selection_intensity = selection_intensity,
    trait_names = trait_names,
    summary = summary_df,
    trait_gains = trait_gains
  )

  class(result) <- c("gw_lmsi", "marker_index", "list")
  result
}


#' @export
print.lmsi <- function(x, ...) {
  cat("\n================================================================\n")
  cat("LINEAR MARKER SELECTION INDEX (LMSI)\n")
  cat("Lande & Thompson (1990) - Chapter 4, Section 4.1\n")
  cat("================================================================\n\n")

  cat("Selection intensity (k):", x$selection_intensity, "\n")
  cat("Number of traits:       ", length(x$trait_names), "\n\n")

  cat("----------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("----------------------------------------------------------------\n")
  metrics <- attr(x$summary, "metrics")
  cat(sprintf("  Accuracy (rHI):         %.4f\n", metrics$rHI))
  cat(sprintf("  Heritability (hI2):     %.4f\n", metrics$hI2))
  cat(sprintf("  Index SD (sigma_I):     %.4f\n", metrics$sigma_I))
  cat(sprintf("  Response (R):           %.4f\n", metrics$R))
  cat(sprintf("  Genetic Advance (GA):   %.4f\n", metrics$GA))
  if (!is.na(metrics$PRE)) {
    cat(sprintf("  Relative Efficiency:    %.2f%%\n", metrics$PRE))
  }

  cat("\n----------------------------------------------------------------\n")
  cat("COEFFICIENT ANALYSIS: PHENOTYPE vs MARKER SCORE WEIGHTING\n")
  cat("----------------------------------------------------------------\n")
  print(x$coeff_analysis, row.names = FALSE, digits = 4)

  cat("\n----------------------------------------------------------------\n")
  cat("PHENOTYPE COEFFICIENTS (b_y)\n")
  cat("----------------------------------------------------------------\n")
  print(x$phenotype_coeffs, row.names = FALSE)

  cat("\n----------------------------------------------------------------\n")
  cat("MARKER SCORE COEFFICIENTS (b_s)\n")
  cat("----------------------------------------------------------------\n")
  print(x$marker_coeffs, row.names = FALSE)

  cat("\n----------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT\n")
  cat("----------------------------------------------------------------\n")
  gains_summary <- data.frame(
    Trait = x$trait_names,
    Weight = x$phenotype_coeffs$weight,
    Delta_H = round(x$Delta_H, 4),
    stringsAsFactors = FALSE
  )
  print(gains_summary, row.names = FALSE)

  cat("\n")
  invisible(x)
}

#' @export
print.gw_lmsi <- function(x, ...) {
  cat("\n================================================================\n")
  cat("GENOME-WIDE LINEAR MARKER SELECTION INDEX (GW-LMSI)\n")
  cat("Lande & Thompson (1990) - Chapter 4, Section 4.2\n")
  cat("================================================================\n\n")

  cat("Selection intensity (k):", x$selection_intensity, "\n")
  cat("Number of traits:       ", x$n_traits, "\n")
  cat("Number of markers:      ", x$n_markers, "\n")
  cat("Number of genotypes:    ", x$n_genotypes, "\n")

  if (x$high_dimensional) {
    cat("Matrix status:           HIGH-DIMENSIONAL (n_markers > n_genotypes)\n")
  } else {
    cat("Matrix status:           Standard (n_markers <= n_genotypes)\n")
  }

  if (!is.na(x$condition_number)) {
    if (x$condition_number > 1e10) {
      cat("Condition number:        ILL-CONDITIONED (", sprintf("%.2e", x$condition_number), ")\n")
    } else if (x$condition_number > 1e6) {
      cat("Condition number:        MODERATE (", sprintf("%.2e", x$condition_number), ")\n")
    } else {
      cat("Condition number:        GOOD (", sprintf("%.2e", x$condition_number), ")\n")
    }
  }

  if (x$ridge_applied) {
    cat("Ridge regularization:    APPLIED (lambda =", x$lambda, ")\n")
  } else {
    cat("Ridge regularization:    NONE (lambda = 0)\n")
  }
  cat("\n")

  cat("----------------------------------------------------------------\n")
  cat("INDEX METRICS\n")
  cat("----------------------------------------------------------------\n")
  print(x$summary, row.names = FALSE)

  cat("\n----------------------------------------------------------------\n")
  cat("EXPECTED GENETIC GAINS PER TRAIT\n")
  cat("----------------------------------------------------------------\n")
  print(x$trait_gains, row.names = FALSE)

  cat("\n")
  invisible(x)
}
