#' Genomic Variance-Covariance Functions
#'
#' @description
#' Functions for computing genomic variance-covariance matrices from GEBVs
#' and combined phenomic-genomic matrices for genomic selection indices.
#'
#' @name genomic-varcov
#' @keywords internal
NULL


#' Genomic Variance-Covariance Matrix (Γ)
#'
#' @description
#' Computes genomic variance-covariance matrix (Γ or Gamma) from a matrix of
#' Genomic Estimated Breeding Values (GEBVs).
#'
#' γ (gamma) represents GEBV vectors obtained from genomic prediction models
#' (e.g., GBLUP, rrBLUP, Genomic BLUP). This function computes Var(γ) = Γ.
#'
#' @param gebv_mat Matrix of GEBVs (n_genotypes x n_traits)
#' @param method Character string specifying correlation method: "pearson" (default),
#'   "kendall", or "spearman"
#' @param use Character string specifying how to handle missing values:
#'   "everything" (default), "complete.obs", "pairwise.complete.obs", etc.
#'   See \code{\link[stats]{cov}} for details.
#'
#' @return Symmetric genomic variance-covariance matrix (n_traits x n_traits)
#'
#' @details
#' The genomic variance-covariance matrix Γ captures genetic variation as
#' predicted by molecular markers. It is computed as:
#'

#'
#' where γ_i is the GEBV vector for genotype i and μ_γ is the mean GEBV vector.
#'
#' **Missing Value Handling:**
#' - "complete.obs": Uses only complete observations (recommended)
#' - "pairwise.complete.obs": Uses pairwise-complete observations (may not be PSD)
#' - "everything": Fails if any NA present
#'
#' When using pairwise deletion, the resulting matrix may not be positive
#' semi-definite (PSD), which can cause numerical issues in selection indices.
#'
#' **Applications:**
#' In selection index theory:
#' - Used in LGSI (Linear Genomic Selection Index)
#' - Component of Φ (phenomic-genomic covariance)
#' - Component of A (genetic-genomic covariance)
#'
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Chapters 4 & 8.
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate GEBVs
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- 5
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits),
#'   nrow = n_genotypes, ncol = n_traits
#' )
#' colnames(gebv_mat) <- paste0("Trait", 1:n_traits)
#'
#' # Compute genomic variance-covariance
#' Gamma <- genomic_varcov(gebv_mat)
#' print(Gamma)
#' }
genomic_varcov <- function(gebv_mat, method = "pearson", use = "complete.obs") {
  gebv_mat <- as.matrix(gebv_mat)
  if (!is.numeric(gebv_mat)) {
    stop("gebv_mat must be numeric")
  }

  n_genotypes <- nrow(gebv_mat)


  if (n_genotypes < 2) {
    stop("gebv_mat must have at least 2 observations")
  }

  has_missing <- any(!is.finite(gebv_mat))
  if (has_missing) {
    if (use == "pairwise.complete.obs") {
      warning("Missing values detected in gebv_mat with use='pairwise.complete.obs'. ",
        "The resulting covariance matrix may not be positive semi-definite, ",
        "which can cause issues in selection index calculations. ",
        "Consider using use='complete.obs' or imputing missing values.",
        call. = FALSE
      )
    } else if (use == "everything") {
      stop(
        "Missing values detected in gebv_mat. Cannot compute covariance with use='everything'. ",
        "Use 'complete.obs', 'pairwise.complete.obs', or impute missing values."
      )
    }
  }

  Gamma <- stats::cov(gebv_mat, use = use, method = method)

  Gamma <- (Gamma + t(Gamma)) / 2

  if (has_missing && use == "pairwise.complete.obs") {
    warn_pairwise_psd(Gamma, "Genomic covariance matrix")
  }

  trait_names <- colnames(gebv_mat)
  if (!is.null(trait_names)) {
    dimnames(Gamma) <- list(trait_names, trait_names)
  }

  Gamma
}


#' Phenomic-Genomic Variance-Covariance Matrix (Φ)
#'
#' @description
#' Computes the combined phenomic-genomic variance-covariance matrix (Φ or P_L),
#' which is the block matrix representing the joint distribution of phenotypes
#' and GEBVs.
#'
#' Structure: Φ = [[P, P_yγ], [P_yγ', Γ]]
#'
#' where:
#' - P = Var(y) = phenotypic variance-covariance
#' - Γ = Var(γ) = genomic variance-covariance
#' - P_yγ = Cov(y, γ) = covariance between phenotypes and GEBVs
#'
#' @param phen_mat Matrix of phenotypes (n_genotypes x n_traits).
#'   Optional if P and P_yg are provided.
#' @param gebv_mat Matrix of GEBVs (n_genotypes x n_traits).
#'   Optional if Gamma and P_yg are provided.
#' @param P Phenotypic variance-covariance matrix (n_traits x n_traits).
#'   Optional if phen_mat is provided.
#' @param Gamma Genomic variance-covariance matrix (n_traits x n_traits).
#'   Optional if gebv_mat is provided.
#' @param P_yg Covariance between phenotypes and GEBVs (n_traits x n_traits).
#'   Optional if phen_mat and gebv_mat are provided.
#' @param method Character string specifying correlation method: "pearson" (default),
#'   "kendall", or "spearman"
#' @param use Character string specifying how to handle missing values:
#'   "complete.obs" (default), "pairwise.complete.obs", etc.
#'
#' @return Symmetric block matrix Φ (2*n_traits x 2*n_traits)
#'
#' @details
#' The phenomic-genomic covariance matrix is used in:
#' - GESIM (Genomic Eigen Selection Index Method)
#' - Combined phenotypic + genomic selection indices
#'
#' The matrix is constructed as:
#' \deqn{\Phi = \begin{bmatrix} P & P_{y\gamma} \\ P_{y\gamma}' & \Gamma \end{bmatrix}}
#'
#' where the off-diagonal blocks are transposes, ensuring symmetry.
#'
#' You can provide either:
#' 1. Raw data: phen_mat + gebv_mat (matrices computed internally)
#' 2. Pre-computed matrices: P + Gamma + P_yg
#'
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Chapter 8.
#'
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' n_genotypes <- 100
#' n_traits <- 7
#' phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
#'   nrow = n_genotypes, ncol = n_traits
#' )
#' gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
#'   nrow = n_genotypes, ncol = n_traits
#' )
#'
#' # Compute phenomic-genomic covariance
#' Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)
#' print(dim(Phi)) # Should be 14 x 14 (2 * 7 traits)
#' }
phenomic_genomic_varcov <- function(phen_mat = NULL, gebv_mat = NULL,
                                    P = NULL, Gamma = NULL, P_yg = NULL,
                                    method = "pearson", use = "complete.obs") {
  has_matrices <- !is.null(P) && !is.null(Gamma) && !is.null(P_yg)
  has_data <- !is.null(phen_mat) && !is.null(gebv_mat)

  if (!has_matrices && !has_data) {
    stop("Must provide either (phen_mat, gebv_mat) or (P, Gamma, P_yg)")
  }

  if (has_data) {
    phen_mat <- as.matrix(phen_mat)
    gebv_mat <- as.matrix(gebv_mat)

    if (nrow(phen_mat) != nrow(gebv_mat)) {
      stop("phen_mat and gebv_mat must have the same number of rows")
    }
    if (ncol(phen_mat) != ncol(gebv_mat)) {
      stop("phen_mat and gebv_mat must have the same number of columns (traits)")
    }

    n_traits <- ncol(phen_mat)

    if (is.null(P)) {
      P <- stats::cov(phen_mat, use = use, method = method)
      P <- (P + t(P)) / 2 # Ensure symmetry
    }

    if (is.null(Gamma)) {
      Gamma <- stats::cov(gebv_mat, use = use, method = method)
      Gamma <- (Gamma + t(Gamma)) / 2 # Ensure symmetry
    }

    if (is.null(P_yg)) {
      P_yg <- stats::cov(phen_mat, gebv_mat, use = use, method = method)
    }
  } else {
    P <- as.matrix(P)
    Gamma <- as.matrix(Gamma)
    P_yg <- as.matrix(P_yg)

    if (!is_symmetric(P)) {
      stop("P must be symmetric")
    }
    if (!is_symmetric(Gamma)) {
      stop("Gamma must be symmetric")
    }

    n_traits <- nrow(P)

    if (nrow(Gamma) != n_traits || ncol(Gamma) != n_traits) {
      stop("Gamma must be ", n_traits, " x ", n_traits)
    }
    if (nrow(P_yg) != n_traits || ncol(P_yg) != n_traits) {
      stop("P_yg must be ", n_traits, " x ", n_traits)
    }
  }

  Phi <- rbind(
    cbind(P, P_yg),
    cbind(t(P_yg), Gamma)
  )

  if (!is_symmetric(Phi, tolerance = TOL_EQUAL)) {
    max_asymmetry <- max(abs(Phi - t(Phi)))
    warning(
      "Phi is not symmetric (max asymmetry = ",
      formatC(max_asymmetry, format = "e", digits = 2),
      "). Check input matrices P, Gamma, and P_yg."
    )
  }

  trait_names <- colnames(Phi)[1:n_traits]
  if (!is.null(trait_names)) {
    all_names <- c(paste0(trait_names, "_phen"), paste0(trait_names, "_gebv"))
    dimnames(Phi) <- list(all_names, all_names)
  }

  Phi
}


#' Genetic-Genomic Variance-Covariance Matrix (A)
#'
#' @description
#' Computes the genetic-genomic covariance matrix (A) as defined in Chapter 8
#' (Equation 8.12) for GESIM and related genomic eigen selection indices.
#'
#' Structure: A = [[C, C_g,γ], [C_γ,g, Γ]]  (2t × 2t, square symmetric)
#'
#' where:
#' - C = Var(g) = true genotypic variance-covariance (t × t)
#' - Γ = Var(γ) = genomic variance-covariance (t × t)
#' - C_g,γ = Cov(g, γ) = covariance between true BVs and GEBVs (t × t)
#' - C_γ,g = Cov(γ, g) = transpose of C_g,γ (t × t)
#'
#' @param gmat Genotypic variance-covariance matrix (n_traits x n_traits)
#' @param Gamma Genomic variance-covariance matrix (n_traits x n_traits).
#'   If NULL, assumed equal to gmat (perfect prediction).
#' @param reliability Optional. Reliability of GEBVs (r² = squared correlation
#'   between GEBV and true BV). Can be:
#'   - Single value (applied to all traits)
#'   - Vector of length n_traits (one per trait)
#'   - NULL (default): assumes C_g,γ = Gamma (unbiased GEBVs with reliability = 1)
#' @param C_gebv_g Optional. Direct specification of Cov(γ, g) matrix (t × t).
#'   If provided, overrides reliability parameter.
#' @param square Logical. If TRUE (default), returns (2t × 2t) square matrix
#'   as required for GESIM. If FALSE, returns (2t × t) rectangular form for LMSI.
#'
#' @return Genetic-genomic covariance matrix:
#'   - If square = TRUE: (2t × 2t) symmetric matrix for GESIM/eigen indices
#'   - If square = FALSE: (2t × t) rectangular matrix for LMSI
#'   where t is the number of traits
#'
#' @details
#' The genetic-genomic matrix relates selection on phenotypes + GEBVs to
#' expected genetic gains.
#'
#' **For GESIM (Chapter 8):** Requires the full (2t × 2t) square matrix for
#' the eigenproblem: (Φ^(-1) A - λI)b = 0
#'
#' **For LMSI/CLGSI (Chapter 4):** Can use the rectangular (2t × t) form
#' in the equation: b = P^(-1) G w, where G is (2t × t).
#'
#' When reliability is provided:
#' - C_γg = diag(√r²) %*% gmat (assumes accuracy scales genetic covariance)
#'
#' When reliability is NULL:
#' - C_γg = Gamma (assumes unbiased GEBVs, perfect prediction)
#'
#' @references
#' Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in Modern
#' Plant Breeding. Springer International Publishing. Chapters 4 & 8.
#'
#' @export
#' @examples
#' \dontrun{
#' # Generate example data
#' gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
#'
#' # Simulate genomic covariance
#' Gamma <- gmat * 0.8
#'
#' # For GESIM: Get square (2t × 2t) matrix
#' A_square <- genetic_genomic_varcov(gmat, Gamma, reliability = 0.7)
#' print(dim(A_square)) # Should be 14 x 14 (2t × 2t)
#'
#' # For LMSI: Get rectangular (2t × t) matrix
#' A_rect <- genetic_genomic_varcov(gmat, Gamma, reliability = 0.7, square = FALSE)
#' print(dim(A_rect)) # Should be 14 x 7 (2t × t)
#' }
genetic_genomic_varcov <- function(gmat, Gamma = NULL, reliability = NULL,
                                   C_gebv_g = NULL, square = TRUE) {
  gmat <- as.matrix(gmat)
  n_traits <- nrow(gmat)

  if (!is_symmetric(gmat)) {
    stop("gmat must be symmetric")
  }

  if (is.null(Gamma)) {
    Gamma <- gmat # Assume perfect prediction
    Gamma <- (Gamma + t(Gamma)) / 2
  } else {
    Gamma <- as.matrix(Gamma)
    if (!is_symmetric(Gamma)) {
      stop("Gamma must be symmetric")
    }
    if (nrow(Gamma) != n_traits || ncol(Gamma) != n_traits) {
      stop("Gamma must be ", n_traits, " x ", n_traits)
    }
  }

  if (!is.null(C_gebv_g)) {
    C_gebv_g <- as.matrix(C_gebv_g)
    if (nrow(C_gebv_g) != n_traits || ncol(C_gebv_g) != n_traits) {
      stop("C_gebv_g must be ", n_traits, " x ", n_traits)
    }
  } else if (!is.null(reliability)) {

    if (length(reliability) == 1) {
      r_squared_vec <- rep(reliability, n_traits)
    } else if (length(reliability) == n_traits) {
      r_squared_vec <- reliability
    } else {
      stop("reliability must be a single value or vector of length ", n_traits)
    }

    if (any(r_squared_vec < 0) || any(r_squared_vec > 1)) {
      stop("reliability values must be between 0 and 1")
    }

    accuracy_vec <- sqrt(r_squared_vec)

    C_gebv_g <- sweep(gmat, 1, accuracy_vec, "*")
  } else {
    C_gebv_g <- Gamma
  }

  if (square) {
    A <- rbind(
      cbind(gmat, C_gebv_g),
      cbind(t(C_gebv_g), Gamma)
    )

    if (!is_symmetric(A, tolerance = TOL_EQUAL)) {
      max_asymmetry <- max(abs(A - t(A)))
      warning(
        "A is not symmetric (max asymmetry = ",
        formatC(max_asymmetry, format = "e", digits = 2),
        "). Check input matrices."
      )
    }

    trait_names <- colnames(gmat)
    if (!is.null(trait_names)) {
      all_names <- c(paste0(trait_names, "_phen"), paste0(trait_names, "_gebv"))
      dimnames(A) <- list(all_names, all_names)
    }
  } else {
    A <- rbind(gmat, C_gebv_g)

    trait_names <- colnames(gmat)
    if (!is.null(trait_names)) {
      row_names <- c(paste0(trait_names, "_phen"), paste0(trait_names, "_gebv"))
      dimnames(A) <- list(row_names, trait_names)
    }
  }

  A
}
