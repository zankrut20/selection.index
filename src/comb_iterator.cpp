#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

//' C++ Iterator for Selection Index Combinations
//' 
//' Vectorized computation of selection indices for all trait combinations.
//' Replaces R loop with single-pass compiled code for massive speedup.
//' 
//' @param pmat Phenotypic variance-covariance matrix (n x n)
//' @param gmat Genotypic variance-covariance matrix (n x n)
//' @param wmat Weight vector (n x 1) or matrix (n x wcol)
//' @param comb_matrix Matrix of trait combinations (ncomb x ncomb_total)
//'   Each column is one combination (0-indexed trait indices)
//' @param wcol Weight column to use (0-indexed)
//' @param const_factor Constant factor (default: 2.063)
//' @param PRE_constant Constant for PRE calculation (default: 100 or 100/GAY)
//' 
//' @return List with:
//'   - IDs: Character vector of combination IDs
//'   - b_matrix: Matrix of b coefficients (ncomb_total x max_traits)
//'   - GAs: Numeric vector of genetic advances
//'   - PREs: Numeric vector of percent relative efficiencies
//' 
//' @keywords internal
// [[Rcpp::export]]
List cpp_comb_iterator(const Eigen::Map<Eigen::MatrixXd>& pmat,
                       const Eigen::Map<Eigen::MatrixXd>& gmat,
                       const Eigen::Map<Eigen::MatrixXd>& wmat,
                       const IntegerMatrix& comb_matrix,
                       int wcol = 0,
                       double const_factor = 2.063,
                       double PRE_constant = 100.0) {
  
  int ncomb = comb_matrix.nrow();          // Number of traits per combination
  int ncomb_total = comb_matrix.ncol();    // Total number of combinations
  
  // Pre-allocate result storage
  std::vector<std::string> IDs(ncomb_total);
  std::vector<std::vector<double>> b_list(ncomb_total);
  NumericVector GAs(ncomb_total);
  NumericVector PREs(ncomb_total);
  
  // Main loop - iterate through all combinations
  for (int combo_idx = 0; combo_idx < ncomb_total; ++combo_idx) {
    
    // Extract trait indices for this combination (convert from 1-indexed R to 0-indexed C++)
    std::vector<int> idx(ncomb);
    std::ostringstream id_stream;
    for (int i = 0; i < ncomb; ++i) {
      idx[i] = comb_matrix(i, combo_idx) - 1;  // Convert to 0-indexed
      if (i > 0) id_stream << ", ";
      id_stream << (idx[i] + 1);  // Store as 1-indexed for display
    }
    IDs[combo_idx] = id_stream.str();
    
    // Extract submatrices for this combination
    MatrixXd p_sub(ncomb, ncomb);
    MatrixXd g_sub(ncomb, ncomb);
    VectorXd w_sub(ncomb);
    
    for (int i = 0; i < ncomb; ++i) {
      w_sub(i) = wmat(idx[i], wcol);
      for (int j = 0; j < ncomb; ++j) {
        p_sub(i, j) = pmat(idx[i], idx[j]);
        g_sub(i, j) = gmat(idx[i], idx[j]);
      }
    }
    
    // Solve: bmat = solve(p_sub, g_sub %*% w_sub)
    // Using Eigen's LDLT solver for symmetric positive definite matrices
    VectorXd gw = g_sub * w_sub;
    VectorXd bmat = p_sub.ldlt().solve(gw);
    
    // Calculate GA: numerator = const_factor * t(bmat) %*% g_sub %*% w_sub
    //               denominator = sqrt(t(bmat) %*% p_sub %*% bmat)
    double numerator = const_factor * bmat.dot(gw);
    double denominator = std::sqrt(bmat.dot(p_sub * bmat));
    double G = numerator / denominator;
    
    // Calculate PRE
    double PRE = G * PRE_constant;
    
    // Store results
    b_list[combo_idx].resize(ncomb);
    for (int i = 0; i < ncomb; ++i) {
      b_list[combo_idx][i] = std::round(bmat(i) * 10000.0) / 10000.0;  // Round to 4 decimals
    }
    GAs[combo_idx] = std::round(G * 10000.0) / 10000.0;
    PREs[combo_idx] = std::round(PRE * 10000.0) / 10000.0;
  }
  
  // Find maximum b vector length for matrix construction
  int max_b_cols = 0;
  for (int i = 0; i < ncomb_total; ++i) {
    if ((int)b_list[i].size() > max_b_cols) {
      max_b_cols = b_list[i].size();
    }
  }
  
  // Construct b_matrix (fill with NA for shorter vectors)
  NumericMatrix b_matrix(ncomb_total, max_b_cols);
  std::fill(b_matrix.begin(), b_matrix.end(), NA_REAL);
  
  for (int i = 0; i < ncomb_total; ++i) {
    for (int j = 0; j < (int)b_list[i].size(); ++j) {
      b_matrix(i, j) = b_list[i][j];
    }
  }
  
  // Set column names for b_matrix
  CharacterVector b_colnames(max_b_cols);
  for (int i = 0; i < max_b_cols; ++i) {
    b_colnames[i] = "b." + std::to_string(i + 1);
  }
  colnames(b_matrix) = b_colnames;
  
  return List::create(
    Named("IDs") = IDs,
    Named("b_matrix") = b_matrix,
    Named("GAs") = GAs,
    Named("PREs") = PREs
  );
}
