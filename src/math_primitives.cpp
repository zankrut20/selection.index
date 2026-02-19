#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

//' Generic C++ Math Primitives for Experimental Design Statistics
//' @name cpp_math_primitives
//'
//' @description
//' Generic mathematical operations optimized with C++/Eigen.
//' No design-specific logic - purely mathematical primitives that can be
//' orchestrated by R code to implement any experimental design.
//'
//' This architecture allows:
//' - Easy addition of new experimental designs (in R only)
//' - C++ speed for heavy computation
//' - Single source of truth (design_stats.R)
//' - Better maintainability and testability

//' Compute Grouped Sums for Matrix Columns
//'
//' @description
//' Efficiently computes grouped sums for all columns of a matrix.
//' Equivalent to rowsum() in R but optimized for multiple columns.
//'
//' @param data_mat Numeric matrix (n_obs x n_traits)
//' @param group_idx Integer vector of group indices (1-based, converted to 0-based internally)
//'
//' @return Matrix of grouped sums (n_groups x n_traits)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_grouped_sums(
    const Eigen::Map<Eigen::MatrixXd>& data_mat,
    const Eigen::Map<Eigen::VectorXi>& group_idx
) {
  const int n_obs = data_mat.rows();
  const int n_traits = data_mat.cols();
  
  // Convert 1-based R indices to 0-based C++
  Eigen::VectorXi groups = group_idx.array() - 1;
  
  // Count unique groups
  int n_groups = groups.maxCoeff() + 1;
  
  // Pre-allocate result matrix
  MatrixXd result(n_groups, n_traits);
  result.setZero();
  
  // Accumulate sums for each group
  for (int i = 0; i < n_obs; ++i) {
    result.row(groups(i)) += data_mat.row(i);
  }
  
  return result;
}

//' Compute Multiple Grouped Sums at Once
//'
//' @description
//' Computes grouped sums for multiple grouping variables simultaneously.
//' More efficient than calling cpp_grouped_sums multiple times.
//'
//' @param data_mat Numeric matrix (n_obs x n_traits)
//' @param group_indices List of integer vectors, each representing a grouping variable
//'
//' @return List of matrices, one for each grouping variable
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List cpp_multi_grouped_sums(
    const Eigen::Map<Eigen::MatrixXd>& data_mat,
    const List& group_indices
) {
  const int n_obs = data_mat.rows();
  const int n_traits = data_mat.cols();
  const int n_groupings = group_indices.size();
  
  List result(n_groupings);
  
  for (int g = 0; g < n_groupings; ++g) {
    Eigen::VectorXi groups = Rcpp::as<Eigen::Map<Eigen::VectorXi>>(group_indices[g]).array() - 1;
    int n_groups = groups.maxCoeff() + 1;
    
    MatrixXd group_sums(n_groups, n_traits);
    group_sums.setZero();
    
    for (int i = 0; i < n_obs; ++i) {
      group_sums.row(groups(i)) += data_mat.row(i);
    }
    
    result[g] = group_sums;
  }
  
  return result;
}

//' Compute Sum of Products Between Grouped Sums
//'
//' @description
//' Efficiently computes sum of products for grouped sum vectors.
//' Equivalent to crossprod(sums1, sums2) in R.
//'
//' @param sums1 Matrix of grouped sums (n_groups x n_traits)
//' @param sums2 Matrix of grouped sums (n_groups x n_traits)
//' @param divisor Scalar to divide sums by (e.g., n_replications)
//'
//' @return Matrix of sum of products (n_traits x n_traits)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_crossprod_divided(
    const Eigen::Map<Eigen::MatrixXd>& sums1,
    const Eigen::Map<Eigen::MatrixXd>& sums2,
    double divisor
) {
  // Compute (t(sums1) %*% sums2) / divisor
  return (sums1.transpose() * sums2) / divisor;
}

//' Compute Correction Factor Matrix
//'
//' @description
//' Computes correction factor for all trait pairs.
//' CF[i,j] = (sum_i * sum_j) / n_obs
//'
//' @param data_mat Numeric matrix (n_obs x n_traits)
//'
//' @return Matrix of correction factors (n_traits x n_traits)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_correction_factor_matrix(
    const Eigen::Map<Eigen::MatrixXd>& data_mat
) {
  const int n_obs = data_mat.rows();
  VectorXd grand_totals = data_mat.colwise().sum();
  
  // Outer product: grand_totals * grand_totals^T / n_obs
  return (grand_totals * grand_totals.transpose()) / n_obs;
}

//' Compute Grand Means
//'
//' @description
//' Computes mean for each trait (column).
//'
//' @param data_mat Numeric matrix (n_obs x n_traits)
//'
//' @return Vector of grand means (n_traits)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd cpp_grand_means(
    const Eigen::Map<Eigen::MatrixXd>& data_mat
) {
  return data_mat.colwise().mean();
}

//' Compute Trait-wise Min/Max
//'
//' @description
//' Computes minimum and maximum for each trait.
//'
//' @param data_mat Numeric matrix (n_obs x n_traits)
//'
//' @return List with 'min' and 'max' vectors
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List cpp_trait_minmax(
    const Eigen::Map<Eigen::MatrixXd>& data_mat
) {
  const int n_traits = data_mat.cols();
  VectorXd mins(n_traits);
  VectorXd maxs(n_traits);
  
  for (int i = 0; i < n_traits; ++i) {
    mins(i) = data_mat.col(i).minCoeff();
    maxs(i) = data_mat.col(i).maxCoeff();
  }
  
  return List::create(
    Named("min") = mins,
    Named("max") = maxs
  );
}

//' Compute Genotype Means Matrix
//'
//' @description
//' Efficiently computes means for each genotype across all traits.
//' Equivalent to rowsum(data, genotypes) / counts but optimized.
//'
//' @param data_mat Numeric matrix (n_obs x n_traits)
//' @param gen_idx Integer vector of genotype indices (1-based)
//'
//' @return Matrix of genotype means (n_genotypes x n_traits)
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_genotype_means(
    const Eigen::Map<Eigen::MatrixXd>& data_mat,
    const Eigen::Map<Eigen::VectorXi>& gen_idx
) {
  const int n_obs = data_mat.rows();
  const int n_traits = data_mat.cols();
  
  // Convert to 0-based
  Eigen::VectorXi groups = gen_idx.array() - 1;
  int n_groups = groups.maxCoeff() + 1;
  
  // Compute sums and counts
  MatrixXd sums(n_groups, n_traits);
  VectorXd counts = VectorXd::Zero(n_groups);
  sums.setZero();
  
  for (int i = 0; i < n_obs; ++i) {
    sums.row(groups(i)) += data_mat.row(i);
    counts(groups(i)) += 1.0;
  }
  
  // Divide each row by its count
  for (int g = 0; g < n_groups; ++g) {
    if (counts(g) > 0) {
      sums.row(g) /= counts(g);
    }
  }
  
  return sums;
}

//' Extract Symmetric Submatrix
//'
//' @description
//' Extracts a symmetric submatrix given row/column indices.
//' Used for selecting trait combinations from full covariance matrices.
//'
//' @param mat Numeric matrix (symmetric)
//' @param indices Integer vector of indices (1-based, converted to 0-based internally)
//'
//' @return Symmetric submatrix
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_extract_submatrix(
  const Eigen::Map<Eigen::MatrixXd>& mat,
  const IntegerVector& indices
) {
  int n = indices.size();
  Eigen::MatrixXd result(n, n);
  
  for (int i = 0; i < n; ++i) {
    int row = indices[i] - 1;  // Convert 1-based to 0-based
    for (int j = 0; j < n; ++j) {
      int col = indices[j] - 1;
      result(i, j) = mat(row, col);
    }
  }
  
  return result;
}

//' Extract Vector Elements
//'
//' @description
//' Extracts specific rows from a column of a matrix.
//' Used for extracting weight vectors for trait combinations.
//'
//' @param mat Numeric matrix
//' @param row_indices Integer vector of row indices (1-based)
//' @param col_index Column index (0-based)
//'
//' @return Vector of extracted elements
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd cpp_extract_vector(
  const Eigen::Map<Eigen::MatrixXd>& mat,
  const IntegerVector& row_indices,
  int col_index
) {
  int n = row_indices.size();
  Eigen::VectorXd result(n);
  
  for (int i = 0; i < n; ++i) {
    result(i) = mat(row_indices[i] - 1, col_index);  // Convert 1-based to 0-based
  }
  
  return result;
}

//' Solve Symmetric Linear System
//'
//' @description
//' Solves Ax = b for symmetric positive definite matrix A using LDLT decomposition.
//' More efficient than general solve() for symmetric matrices.
//'
//' @param A Symmetric positive definite matrix
//' @param b Right-hand side vector
//'
//' @return Solution vector x
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd cpp_symmetric_solve(
  const Eigen::Map<Eigen::MatrixXd>& A,
  const Eigen::Map<Eigen::VectorXd>& b
) {
  // Use LDLT decomposition for symmetric matrices
  return A.ldlt().solve(b);
}

//' Quadratic Form: x' A y
//'
//' @description
//' Computes the quadratic form x' A y efficiently.
//' Equivalent to t(x) %*% A %*% y in R but optimized.
//'
//' @param x First vector
//' @param A Matrix
//' @param y Second vector
//'
//' @return Scalar result of x' A y
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double cpp_quadratic_form(
  const Eigen::Map<Eigen::VectorXd>& x,
  const Eigen::Map<Eigen::MatrixXd>& A,
  const Eigen::Map<Eigen::VectorXd>& y
) {
  return x.dot(A * y);
}

//' Symmetric Quadratic Form: x' A x
//'
//' @description
//' Computes the symmetric quadratic form x' A x efficiently.
//' Equivalent to t(x) %*% A %*% x in R but optimized.
//'
//' @param x Vector
//' @param A Symmetric matrix
//'
//' @return Scalar result of x' A x
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double cpp_quadratic_form_sym(
  const Eigen::Map<Eigen::VectorXd>& x,
  const Eigen::Map<Eigen::MatrixXd>& A
) {
  return x.dot(A * x);
}

//' @title Correction Factor Matrix
//'
//' @description
//' Computes the correction factor matrix for ANOVA calculations.
//' CF(i,j) = (sum_i * sum_j) / n
//'
//' @param total_sums Vector of column sums
//' @param n_obs Number of observations
//'
//' @return Symmetric correction factor matrix
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_correction_factor(
  const Eigen::Map<Eigen::VectorXd>& total_sums,
  int n_obs
) {
  int n_traits = total_sums.size();
  Eigen::MatrixXd CF(n_traits, n_traits);
  
  for (int i = 0; i < n_traits; ++i) {
    for (int j = i; j < n_traits; ++j) {
      CF(i, j) = (total_sums[i] * total_sums[j]) / n_obs;
      if (i != j) CF(j, i) = CF(i, j);
    }
  }
  
  return CF;
}

//' @title Total Sum of Products
//'
//' @description
//' Computes the total sum of products matrix corrected for the mean.
//' TSP(i,j) = sum(x_i * x_j) - CF(i,j)
//'
//' @param data_mat Data matrix (n_obs x n_traits)
//' @param CF Correction factor matrix
//'
//' @return Symmetric sum of products matrix
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_total_sum_of_products(
  const Eigen::Map<Eigen::MatrixXd>& data_mat,
  const Eigen::Map<Eigen::MatrixXd>& CF
) {
  int n_traits = data_mat.cols();
  Eigen::MatrixXd TSP(n_traits, n_traits);
  
  for (int i = 0; i < n_traits; ++i) {
    for (int j = i; j < n_traits; ++j) {
      TSP(i, j) = data_mat.col(i).dot(data_mat.col(j)) - CF(i, j);
      if (i != j) TSP(j, i) = TSP(i, j);
    }
  }
  
  return TSP;
}

//' @title Grouped Sum of Products
//'
//' @description
//' Computes the sum of products for grouped data.
//' GSP(i,j) = sum_g [(sum_i_g * sum_j_g) / n_g] - CF(i,j)
//'
//' @param group_sums Matrix of group sums (n_groups x n_traits)
//' @param group_counts Vector of group sizes
//' @param CF Correction factor matrix
//'
//' @return Symmetric sum of products matrix
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_grouped_sum_of_products(
  const Eigen::Map<Eigen::MatrixXd>& group_sums,
  const Eigen::Map<Eigen::VectorXi>& group_counts,
  const Eigen::Map<Eigen::MatrixXd>& CF
) {
  int n_groups = group_sums.rows();
  int n_traits = group_sums.cols();
  Eigen::MatrixXd GSP = Eigen::MatrixXd::Zero(n_traits, n_traits);
  
  for (int g = 0; g < n_groups; ++g) {
    for (int i = 0; i < n_traits; ++i) {
      for (int j = i; j < n_traits; ++j) {
        GSP(i, j) += (group_sums(g, i) * group_sums(g, j)) / group_counts[g];
        if (i != j) GSP(j, i) = GSP(i, j);
      }
    }
  }
  
  // Subtract correction factor
  GSP -= CF;
  
  return GSP;
}

//' @title Mean Squares from Sum of Products
//'
//' @description
//' Divides sum of products matrix by degrees of freedom.
//' MS = SP / df
//'
//' @param sum_of_products Sum of products matrix
//' @param df Degrees of freedom
//'
//' @return Mean squares matrix
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cpp_mean_squares(
  const Eigen::Map<Eigen::MatrixXd>& sum_of_products,
  int df
) {
  return sum_of_products / df;
}
