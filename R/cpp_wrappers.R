#' C++ Function Wrappers with Uniform Validation
#'
#' @description
#' Thin R wrapper functions that call cpp_* functions with uniform input validation,
#' type checking, and error handling. Centralizes validation logic and provides
#' consistent error messages across the package.
#'
#' @name cpp_wrappers
#' @keywords internal
NULL


#' Grouped Sums Wrapper
#'
#' @description
#' Wrapper for cpp_grouped_sums() with validation.
#' Computes grouped sums for all matrix columns efficiently.
#'
#' @param data_mat Numeric matrix (n_obs x n_traits)
#' @param group_idx Integer vector of group indices
#' @param check_na Logical; check for NA values (default: TRUE)
#'
#' @return Matrix of grouped sums (n_groups x n_traits)
#'
#' @keywords internal
#' @noRd
grouped_sums <- function(data_mat, group_idx, check_na = TRUE) {
  
  # Validate data matrix
  if (!is.matrix(data_mat)) {
    data_mat <- as.matrix(data_mat)
  }
  if (!is.numeric(data_mat)) {
    stop("data_mat must be numeric")
  }
  storage.mode(data_mat) <- "numeric"
  
  # Check for NA values if requested
  if (check_na && anyNA(data_mat)) {
    stop("data_mat contains NA values")
  }
  
  # Validate group indices
  if (!is.vector(group_idx)) {
    stop("group_idx must be a vector")
  }
  if (length(group_idx) != nrow(data_mat)) {
    stop("Length of group_idx (", length(group_idx), ") must match ",
         "number of rows in data_mat (", nrow(data_mat), ")")
  }
  if (anyNA(group_idx)) {
    stop("group_idx contains NA values")
  }
  
  # Convert to integer if needed
  if (!is.integer(group_idx)) {
    group_idx <- as.integer(group_idx)
  }
  
  # Call C++ function
  cpp_grouped_sums(data_mat, group_idx)
}


#' Correction Factor Wrapper
#'
#' @description
#' Wrapper for cpp_correction_factor() with validation.
#' Computes correction factor matrix for ANOVA.
#'
#' @param total_sums Numeric vector of column sums
#' @param n_obs Integer; number of observations
#'
#' @return Symmetric correction factor matrix
#'
#' @keywords internal
#' @noRd
correction_factor <- function(total_sums, n_obs) {
  
  # Validate total sums
  if (!is.numeric(total_sums)) {
    stop("total_sums must be numeric")
  }
  storage.mode(total_sums) <- "numeric"
  
  if (anyNA(total_sums)) {
    stop("total_sums contains NA values")
  }
  
  # Validate n_obs
  if (!is.numeric(n_obs) || length(n_obs) != 1) {
    stop("n_obs must be a single numeric value")
  }
  if (n_obs < 1) {
    stop("n_obs must be positive (got ", n_obs, ")")
  }
  n_obs <- as.integer(n_obs)
  
  # Call C++ function
  cpp_correction_factor(total_sums, n_obs)
}


#' Total Sum of Products Wrapper
#'
#' @description
#' Wrapper for cpp_total_sum_of_products() with validation.
#' Computes total sum of products matrix corrected for mean.
#'
#' @param data_mat Numeric matrix (n_obs x n_traits)
#' @param CF Correction factor matrix
#'
#' @return Symmetric sum of products matrix
#'
#' @keywords internal
#' @noRd
total_sum_of_products <- function(data_mat, CF) {
  
  # Validate data matrix
  if (!is.matrix(data_mat)) {
    data_mat <- as.matrix(data_mat)
  }
  if (!is.numeric(data_mat)) {
    stop("data_mat must be numeric")
  }
  storage.mode(data_mat) <- "numeric"
  
  # Validate correction factor
  if (!is.matrix(CF)) {
    stop("CF must be a matrix")
  }
  if (!is.numeric(CF)) {
    stop("CF must be numeric")
  }
  
  # Check dimensions match
  n_traits <- ncol(data_mat)
  if (nrow(CF) != n_traits || ncol(CF) != n_traits) {
    stop("CF dimensions (", nrow(CF), "x", ncol(CF), ") must match ",
         "number of traits (", n_traits, ")")
  }
  
  # Call C++ function
  cpp_total_sum_of_products(data_mat, CF)
}


#' Grouped Sum of Products Wrapper
#'
#' @description
#' Wrapper for cpp_grouped_sum_of_products() with validation.
#' Computes sum of products for grouped data.
#'
#' @param group_sums Matrix of group sums (n_groups x n_traits)
#' @param group_counts Integer vector of group sizes
#' @param CF Correction factor matrix
#'
#' @return Symmetric sum of products matrix
#'
#' @keywords internal
#' @noRd
grouped_sum_of_products <- function(group_sums, group_counts, CF) {
  
  # Validate group sums
  if (!is.matrix(group_sums)) {
    group_sums <- as.matrix(group_sums)
  }
  if (!is.numeric(group_sums)) {
    stop("group_sums must be numeric")
  }
  storage.mode(group_sums) <- "numeric"
  
  # Validate group counts
  if (!is.vector(group_counts)) {
    stop("group_counts must be a vector")
  }
  if (length(group_counts) != nrow(group_sums)) {
    stop("Length of group_counts (", length(group_counts), ") must match ",
         "number of groups in group_sums (", nrow(group_sums), ")")
  }
  if (!is.integer(group_counts)) {
    group_counts <- as.integer(group_counts)
  }
  if (any(group_counts < 1)) {
    stop("All group_counts must be positive")
  }
  
  # Validate correction factor
  if (!is.matrix(CF)) {
    stop("CF must be a matrix")
  }
  if (!is.numeric(CF)) {
    stop("CF must be numeric")
  }
  
  # Check dimensions match
  n_traits <- ncol(group_sums)
  if (nrow(CF) != n_traits || ncol(CF) != n_traits) {
    stop("CF dimensions (", nrow(CF), "x", ncol(CF), ") must match ",
         "number of traits (", n_traits, ")")
  }
  
  # Call C++ function
  cpp_grouped_sum_of_products(group_sums, group_counts, CF)
}


#' Mean Squares Wrapper
#'
#' @description
#' Wrapper for cpp_mean_squares() with validation.
#' Divides sum of products matrix by degrees of freedom.
#'
#' @param sum_of_products Sum of products matrix
#' @param df Integer; degrees of freedom
#'
#' @return Mean squares matrix
#'
#' @keywords internal
#' @noRd
mean_squares <- function(sum_of_products, df) {
  
  # Validate sum of products
  if (!is.matrix(sum_of_products)) {
    sum_of_products <- as.matrix(sum_of_products)
  }
  if (!is.numeric(sum_of_products)) {
    stop("sum_of_products must be numeric")
  }
  
  # Validate degrees of freedom
  if (!is.numeric(df) || length(df) != 1) {
    stop("df must be a single numeric value")
  }
  if (df < 1) {
    stop("df must be positive (got ", df, ")")
  }
  df <- as.integer(df)
  
  # Call C++ function
  cpp_mean_squares(sum_of_products, df)
}


#' Genotype Means Wrapper
#'
#' @description
#' Wrapper for cpp_genotype_means() with validation.
#' Efficiently computes means for each genotype across all traits.
#'
#' @param data_mat Numeric matrix (n_obs x n_traits)
#' @param gen_idx Integer vector of genotype indices
#' @param check_na Logical; check for NA values (default: TRUE)
#'
#' @return Matrix of genotype means (n_genotypes x n_traits)
#'
#' @keywords internal
#' @noRd
genotype_means <- function(data_mat, gen_idx, check_na = TRUE) {
  
  # Validate data matrix
  if (!is.matrix(data_mat)) {
    data_mat <- as.matrix(data_mat)
  }
  if (!is.numeric(data_mat)) {
    stop("data_mat must be numeric")
  }
  storage.mode(data_mat) <- "numeric"
  
  # Check for NA values if requested
  if (check_na && anyNA(data_mat)) {
    stop("data_mat contains NA values")
  }
  
  # Validate genotype indices
  if (!is.vector(gen_idx)) {
    stop("gen_idx must be a vector")
  }
  if (length(gen_idx) != nrow(data_mat)) {
    stop("Length of gen_idx (", length(gen_idx), ") must match ",
         "number of rows in data_mat (", nrow(data_mat), ")")
  }
  if (anyNA(gen_idx)) {
    stop("gen_idx contains NA values")
  }
  
  # Convert to integer if needed
  if (!is.integer(gen_idx)) {
    gen_idx <- as.integer(gen_idx)
  }
  
  # Call C++ function
  cpp_genotype_means(data_mat, gen_idx)
}


#' Symmetric Solve Wrapper
#'
#' @description
#' Wrapper for cpp_symmetric_solve() with validation.
#' Solves Ax = b for symmetric positive definite matrix A.
#'
#' @param A Symmetric positive definite matrix
#' @param b Right-hand side vector or matrix
#'
#' @return Solution vector/matrix x
#'
#' @keywords internal
#' @noRd
symmetric_solve <- function(A, b) {
  
  # Validate A
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  if (!is.numeric(A)) {
    stop("A must be numeric")
  }
  
  # Check A is square
  if (nrow(A) != ncol(A)) {
    stop("A must be square (got ", nrow(A), "x", ncol(A), ")")
  }
  
  # Check symmetry
  if (!is_symmetric(A)) {
    warning("A is not symmetric (within tolerance ", TOL_SYM, "). ",
            "Results may be unreliable.")
  }
  
  # Validate b
  if (is.matrix(b)) {
    if (!is.numeric(b)) {
      stop("b must be numeric")
    }
    if (nrow(b) != nrow(A)) {
      stop("Number of rows in b (", nrow(b), ") must match dimension of A (", nrow(A), ")")
    }
  } else {
    if (!is.numeric(b)) {
      stop("b must be numeric")
    }
    if (length(b) != nrow(A)) {
      stop("Length of b (", length(b), ") must match dimension of A (", nrow(A), ")")
    }
  }
  
  # Call C++ function
  cpp_symmetric_solve(A, b)
}


#' Quadratic Form Wrapper
#'
#' @description
#' Wrapper for cpp_quadratic_form() with validation.
#' Computes x' A y efficiently.
#'
#' @param x First vector
#' @param A Matrix
#' @param y Second vector
#'
#' @return Scalar result of x' A y
#'
#' @keywords internal
#' @noRd
quadratic_form <- function(x, A, y) {
  
  # Validate x
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  if (!is.vector(x)) {
    x <- as.vector(x)
  }
  storage.mode(x) <- "numeric"
  
  # Validate A
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  if (!is.numeric(A)) {
    stop("A must be numeric")
  }
  
  # Validate y
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (!is.vector(y)) {
    y <- as.vector(y)
  }
  storage.mode(y) <- "numeric"
  
  # Check dimensions
  if (length(x) != nrow(A)) {
    stop("Length of x (", length(x), ") must match rows of A (", nrow(A), ")")
  }
  if (length(y) != ncol(A)) {
    stop("Length of y (", length(y), ") must match columns of A (", ncol(A), ")")
  }
  
  # Call C++ function
  cpp_quadratic_form(x, A, y)
}


#' Symmetric Quadratic Form Wrapper
#'
#' @description
#' Wrapper for cpp_quadratic_form_sym() with validation.
#' Computes x' A x efficiently.
#'
#' @param x Vector
#' @param A Symmetric matrix
#'
#' @return Scalar result of x' A x
#'
#' @keywords internal
#' @noRd
quadratic_form_sym <- function(x, A) {
  
  # Validate x
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  if (!is.vector(x)) {
    x <- as.vector(x)
  }
  storage.mode(x) <- "numeric"
  
  # Validate A
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  if (!is.numeric(A)) {
    stop("A must be numeric")
  }
  
  # Check A is square
  if (nrow(A) != ncol(A)) {
    stop("A must be square (got ", nrow(A), "x", ncol(A), ")")
  }
  
  # Check dimension compatibility
  if (length(x) != nrow(A)) {
    stop("Length of x (", length(x), ") must match dimension of A (", nrow(A), ")")
  }
  
  # Call C++ function
  cpp_quadratic_form_sym(x, A)
}
