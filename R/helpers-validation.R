#' Validation Helper Functions
#'
#' @description
#' Internal validation functions for experimental design arguments,
#' index vectors, and matrix properties. Centralizes validation logic
#' to eliminate duplication across the package.
#'
#' @name helpers-validation
#' @keywords internal
#' @noRd
NULL


#' Validate Design Arguments
#'
#' @description
#' Validates experimental design type and required parameters.
#' Ensures design-specific requirements are met (e.g., columns for LSD,
#' main_plots for SPD).
#'
#' @param design_type Integer design code (1=RCBD, 2=LSD, 3=SPD) or
#'   character string ("RCBD", "LSD", "SPD")
#' @param col_idx Integer vector of column indices (required for LSD)
#' @param main_idx Integer vector of main plot indices (required for SPD)
#' @param allow_char Logical; if TRUE, allow character design_type and
#'   convert to integer code
#'
#' @return Integer design code (1, 2, or 3) if validation passes
#'
#' @details
#' Stops with informative error if:
#' - design_type is not valid (must be 1, 2, 3 or "RCBD", "LSD", "SPD")
#' - LSD is specified but col_idx is NULL
#' - SPD is specified but main_idx is NULL
#'
#' @keywords internal
#' @noRd
validate_design_args <- function(design_type, col_idx = NULL, main_idx = NULL,
                                 allow_char = FALSE) {
  # Convert character design_type to integer code if allowed
  if (allow_char && is.character(design_type)) {
    design_type <- switch(design_type,
      "RCBD" = DESIGN_RCBD,
      "LSD" = DESIGN_LSD,
      "SPD" = DESIGN_SPD,
      stop("design_type must be 'RCBD', 'LSD', or 'SPD'")
    )
  }

  # Validate design_type is a valid code
  if (!design_type %in% DESIGN_CODES) {
    stop("design_type must be 1 (RCBD), 2 (LSD), or 3 (SPD)")
  }

  # Validate LSD requirements
  if (design_type == DESIGN_LSD && is.null(col_idx)) {
    stop("Latin Square Design (design_type = 2) requires 'col_idx' parameter")
  }

  # Validate SPD requirements
  if (design_type == DESIGN_SPD && is.null(main_idx)) {
    stop("Split Plot Design (design_type = 3) requires 'main_idx' parameter")
  }

  design_type
}


#' Validate Index Vectors
#'
#' @description
#' Validates that index vectors (genotypes, replications, columns, main plots)
#' have the correct length matching the data dimensions.
#'
#' @param n_obs Integer; number of observations (rows in data)
#' @param gen_idx Integer vector of genotype indices
#' @param rep_idx Integer vector of replication indices
#' @param col_idx Integer vector of column indices (optional, for LSD)
#' @param main_idx Integer vector of main plot indices (optional, for SPD)
#' @param data_name Character; name of the data object for error messages
#'
#' @return NULL (invisibly) if validation passes; stops with error otherwise
#'
#' @details
#' Checks that:
#' - All index vectors have length equal to n_obs
#' - Index vectors contain at least 2 unique levels (for ANOVA)
#' - No NA values in index vectors
#'
#' @keywords internal
#' @noRd
validate_indices <- function(n_obs, gen_idx, rep_idx,
                             col_idx = NULL, main_idx = NULL,
                             data_name = "data") {
  # Validate genotype indices
  if (length(gen_idx) != n_obs) {
    stop(
      "Length of 'genotypes' (", length(gen_idx), ") must match ",
      "number of rows in ", data_name, " (", n_obs, ")"
    )
  }

  if (any(is.na(gen_idx))) {
    stop("'genotypes' contains NA values")
  }

  n_gen <- length(unique(gen_idx))
  if (n_gen < 2) {
    stop("'genotypes' must have at least 2 unique levels (found ", n_gen, ")")
  }

  # Validate replication indices
  if (length(rep_idx) != n_obs) {
    stop(
      "Length of 'replications' (", length(rep_idx), ") must match ",
      "number of rows in ", data_name, " (", n_obs, ")"
    )
  }

  if (any(is.na(rep_idx))) {
    stop("'replications' contains NA values")
  }

  n_rep <- length(unique(rep_idx))
  if (n_rep < 2) {
    stop("'replications' must have at least 2 unique levels (found ", n_rep, ")")
  }

  # Validate column indices (LSD)
  if (!is.null(col_idx)) {
    if (length(col_idx) != n_obs) {
      stop(
        "Length of 'columns' (", length(col_idx), ") must match ",
        "number of rows in ", data_name, " (", n_obs, ")"
      )
    }

    if (any(is.na(col_idx))) {
      stop("'columns' contains NA values")
    }

    n_col <- length(unique(col_idx))
    if (n_col < 2) {
      stop("'columns' must have at least 2 unique levels (found ", n_col, ")")
    }
  }

  # Validate main plot indices (SPD)
  if (!is.null(main_idx)) {
    if (length(main_idx) != n_obs) {
      stop(
        "Length of 'main_plots' (", length(main_idx), ") must match ",
        "number of rows in ", data_name, " (", n_obs, ")"
      )
    }

    if (any(is.na(main_idx))) {
      stop("'main_plots' contains NA values")
    }

    n_main <- length(unique(main_idx))
    if (n_main < 2) {
      stop("'main_plots' must have at least 2 unique levels (found ", n_main, ")")
    }
  }

  invisible(NULL)
}


#' Warn About Non-Positive Semi-Definite Matrices (Pairwise)
#'
#' @description
#' Checks if a variance-covariance matrix is positive semi-definite (PSD)
#' by examining eigenvalues. Issues a warning if the matrix is not PSD,
#' which can occur with missing data imputation or numerical issues.
#'
#' @param mat Numeric matrix to check (should be symmetric)
#' @param mat_name Character; name of the matrix for warning messages
#' @param tolerance Numeric; tolerance for considering eigenvalues non-negative
#'   (default: TOL_PSD = 1e-8)
#' @param check_symmetry Logical; if TRUE, also check matrix symmetry
#'
#' @return Logical; TRUE if matrix is PSD (within tolerance), FALSE otherwise
#'   Also issues a warning if not PSD.
#'
#' @details
#' A matrix is positive semi-definite if all eigenvalues are non-negative.
#' This function:
#' 1. Optionally checks symmetry first (required for PSD)
#' 2. Computes eigenvalues
#' 3. Checks if minimum eigenvalue >= -tolerance
#' 4. Warns if matrix is not PSD with diagnostic information
#'
#' Common causes of non-PSD matrices:
#' - Missing value imputation artifacts
#' - Numerical precision issues
#' - Invalid covariance estimates
#'
#' @keywords internal
#' @noRd
warn_pairwise_psd <- function(mat, mat_name = "Matrix",
                              tolerance = TOL_PSD,
                              check_symmetry = TRUE) {
  # Check dimensions
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }

  if (nrow(mat) != ncol(mat)) {
    warning(mat_name, " is not square (", nrow(mat), " x ", ncol(mat), ")")
    return(FALSE)
  }

  # Check symmetry if requested
  if (check_symmetry && !isSymmetric(unname(mat), tol = TOL_SYM)) {
    warning(mat_name, " is not symmetric (within tolerance ", TOL_SYM, ")")
    return(FALSE)
  }

  # Compute eigenvalues (symmetric = TRUE for efficiency)
  eigen_vals <- tryCatch(
    eigen(mat, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) {
      warning(mat_name, " eigenvalue computation failed: ", e$message)
      NULL
    }
  )

  if (is.null(eigen_vals)) {
    return(FALSE)
  }

  # Check for positive semi-definiteness
  min_eigenvalue <- min(eigen_vals)

  if (min_eigenvalue < -tolerance) {
    warning(
      mat_name, " is not positive semi-definite.\n",
      "  Minimum eigenvalue: ", formatC(min_eigenvalue, format = "e", digits = 2), "\n",
      "  Maximum eigenvalue: ", formatC(max(eigen_vals), format = "e", digits = 2), "\n",
      "  This may occur with missing data or numerical precision issues.\n",
      "  Consider using a different imputation method or checking data quality.",
      call. = FALSE
    )
    return(FALSE)
  }

  TRUE
}


#' Check Matrix Symmetry
#'
#' @description
#' Convenience wrapper around isSymmetric with package-standard tolerance.
#'
#' @param mat Numeric matrix to check
#' @param tolerance Numeric; tolerance for symmetry check (default: TOL_SYM)
#'
#' @return Logical; TRUE if matrix is symmetric within tolerance
#'
#' @keywords internal
#' @noRd
is_symmetric <- function(mat, tolerance = TOL_SYM) {
  isSymmetric(unname(mat), tol = tolerance)
}


#' Check if Value is Effectively Zero
#'
#' @description
#' Checks if a numeric value is effectively zero within numerical tolerance.
#'
#' @param x Numeric value or vector
#' @param tolerance Numeric; tolerance for zero check (default: TOL_ZERO)
#'
#' @return Logical; TRUE if |x| < tolerance
#'
#' @keywords internal
#' @noRd
is_zero <- function(x, tolerance = TOL_ZERO) {
  # Handle NA and Inf values
  result <- abs(x) < tolerance
  result[is.na(x) | is.infinite(x)] <- FALSE
  result
}
