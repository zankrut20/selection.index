#' Package Constants
#'
#' @description
#' Central location for package-wide constants including design codes,
#' numerical tolerances, and other configuration values.
#'
#' @name constants
#' @keywords internal
NULL

# ============================================================================
# Design Type Codes
# ============================================================================

#' Design Type Code: Randomized Complete Block Design (RCBD)
#' @keywords internal
DESIGN_RCBD <- 1L

#' Design Type Code: Latin Square Design (LSD)
#' @keywords internal
DESIGN_LSD <- 2L

#' Design Type Code: Split Plot Design (SPD)
#' @keywords internal
DESIGN_SPD <- 3L

#' Vector of all valid design codes
#' @keywords internal
DESIGN_CODES <- c(DESIGN_RCBD, DESIGN_LSD, DESIGN_SPD)

#' Named vector mapping design names to codes
#' @keywords internal
DESIGN_NAMES <- c("RCBD" = DESIGN_RCBD, "LSD" = DESIGN_LSD, "SPD" = DESIGN_SPD)


# ============================================================================
# Numerical Tolerances
# ============================================================================

#' Tolerance for symmetry checks (e.g., isSymmetric)
#' Used when checking if matrices are symmetric within numerical precision
#' @keywords internal
TOL_SYM <- 1e-8

#' Tolerance for positive semi-definite checks
#' Eigenvalues greater than -TOL_PSD are considered non-negative
#' @keywords internal
TOL_PSD <- 1e-8

#' Tolerance for zero checks (near-zero values)
#' Values smaller than this in absolute value are considered zero
#' @keywords internal
TOL_ZERO <- 1e-10

#' Tolerance for convergence in iterative algorithms
#' Default convergence criterion for missing value imputation
#' @keywords internal
TOL_CONV <- 1e-6

#' Tolerance for general floating point comparisons
#' @keywords internal
TOL_EQUAL <- 1e-6


# ============================================================================
# Variance/Covariance Type Codes
# ============================================================================

#' Covariance Type Code: Genotypic variance-covariance
#' @keywords internal
COV_GENOTYPIC <- 1L

#' Covariance Type Code: Phenotypic variance-covariance
#' @keywords internal
COV_PHENOTYPIC <- 2L
