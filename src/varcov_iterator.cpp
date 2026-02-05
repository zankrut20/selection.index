#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

//' C++ Variance-Covariance Matrix Iterator
//'
//' @description
//' Optimized C++ implementation for computing variance-covariance matrices.
//' Replaces nested R loops in gen_varcov() and phen_varcov() with vectorized
//' operations. Processes all trait pairs in a single call using Eigen for
//' efficient matrix operations.
//'
//' @param data_mat Numeric matrix of trait observations (n_obs x n_traits)
//' @param gen_idx Integer vector of genotype indices (1-based converted to 0-based internally)
//' @param rep_idx Integer vector of replication indices (1-based converted to 0-based internally)
//' @param col_idx Integer vector of column indices for LSD (NULL for RCBD/SPD)
//' @param main_idx Integer vector of main plot indices for SPD (NULL for RCBD/LSD)
//' @param design_type Integer: 1=RCBD, 2=LSD, 3=SPD
//' @param cov_type Integer: 1=genotypic, 2=phenotypic
//'
//' @return Symmetric variance-covariance matrix (n_traits x n_traits)
//'
//' @details
//' Algorithm optimizations:
//' - Single pass through data for all trait pairs (eliminates R loop overhead)
//' - Uses Eigen's efficient matrix operations (crossprod, grouped sums)
//' - Pre-allocates all result matrices (no memory reallocation)
//' - Computes only upper triangle, then mirrors to lower (symmetric matrix)
//'
//' RCBD: Cov_G = (GMP - EMP) / r, Cov_P = Cov_G + EMP
//' LSD:  Cov_G = (GMP - EMP) / t, Cov_P = Cov_G + EMP  
//' SPD:  Cov_G = (GMP - EMP) / (r*a), Cov_P = Cov_G + EMP
//'
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd cpp_varcov_iterator(
    const Eigen::Map<Eigen::MatrixXd>& data_mat,
    const Eigen::Map<Eigen::VectorXi>& gen_idx,
    const Eigen::Map<Eigen::VectorXi>& rep_idx,
    Nullable<Eigen::Map<Eigen::VectorXi>> col_idx = R_NilValue,
    Nullable<Eigen::Map<Eigen::VectorXi>> main_idx = R_NilValue,
    int design_type = 1,  // 1=RCBD, 2=LSD, 3=SPD
    int cov_type = 1      // 1=genotypic, 2=phenotypic
) {
  
  const int n_obs = data_mat.rows();
  const int n_traits = data_mat.cols();
  
  // Convert 1-based R indices to 0-based C++ indices
  Eigen::VectorXi gen = gen_idx.array() - 1;
  Eigen::VectorXi rep = rep_idx.array() - 1;
  
  // Count unique levels
  int n_gen = gen.maxCoeff() + 1;
  int n_rep = rep.maxCoeff() + 1;
  int n_col = 0;
  int n_main = 0;
  
  Eigen::VectorXi col;
  Eigen::VectorXi main;
  
  if (design_type == 2 && col_idx.isNotNull()) {
    col = Rcpp::as<Eigen::Map<Eigen::VectorXi>>(col_idx).array() - 1;
    n_col = col.maxCoeff() + 1;
  }
  
  if (design_type == 3 && main_idx.isNotNull()) {
    main = Rcpp::as<Eigen::Map<Eigen::VectorXi>>(main_idx).array() - 1;
    n_main = main.maxCoeff() + 1;
  }
  
  // Pre-allocate result matrix (symmetric)
  MatrixXd result = MatrixXd::Zero(n_traits, n_traits);
  
  // Pre-allocate grouped sum matrices for all traits at once
  MatrixXd gen_sums(n_gen, n_traits);
  MatrixXd rep_sums(n_rep, n_traits);
  MatrixXd col_sums, main_sums;
  
  if (design_type == 2) {
    col_sums.resize(n_col, n_traits);
    col_sums.setZero();
  }
  if (design_type == 3) {
    main_sums.resize(n_main, n_traits);
    main_sums.setZero();
  }
  
  // OPTIMIZATION: Compute all grouped sums in one pass
  gen_sums.setZero();
  rep_sums.setZero();
  
  for (int i = 0; i < n_obs; ++i) {
    gen_sums.row(gen(i)) += data_mat.row(i);
    rep_sums.row(rep(i)) += data_mat.row(i);
    
    if (design_type == 2) {
      col_sums.row(col(i)) += data_mat.row(i);
    }
    if (design_type == 3) {
      main_sums.row(main(i)) += data_mat.row(i);
    }
  }
  
  // Grand totals for all traits
  VectorXd grand_totals = data_mat.colwise().sum();
  
  // Process all trait pairs (upper triangle + diagonal)
  for (int i = 0; i < n_traits; ++i) {
    for (int j = i; j < n_traits; ++j) {
      
      // Get trait vectors
      VectorXd trait1 = data_mat.col(i);
      VectorXd trait2 = data_mat.col(j);
      
      // Correction Factor
      double CF = (grand_totals(i) * grand_totals(j)) / n_obs;
      
      // Total Sum of Products
      double TSP = trait1.dot(trait2) - CF;
      
      // Genotype Sum of Products
      VectorXd gen_sum1 = gen_sums.col(i);
      VectorXd gen_sum2 = gen_sums.col(j);
      double GSP = gen_sum1.dot(gen_sum2) / n_rep - CF;
      
      // Replication Sum of Products
      VectorXd rep_sum1 = rep_sums.col(i);
      VectorXd rep_sum2 = rep_sums.col(j);
      double RSP = rep_sum1.dot(rep_sum2) / n_gen - CF;
      
      double GMP, EMP, covariance;
      
      if (design_type == 1) {
        // ========== RCBD ==========
        double ESP = TSP - GSP - RSP;
        int DFE = (n_gen - 1) * (n_rep - 1);
        
        GMP = GSP / (n_gen - 1);
        EMP = ESP / DFE;
        
        if (cov_type == 1) {
          // Genotypic covariance = (GMP - EMP) / r
          covariance = (GMP - EMP) / n_rep;
        } else {
          // Phenotypic covariance = Genotypic + EMP
          covariance = (GMP - EMP) / n_rep + EMP;
        }
        
      } else if (design_type == 2) {
        // ========== Latin Square Design ==========
        VectorXd col_sum1 = col_sums.col(i);
        VectorXd col_sum2 = col_sums.col(j);
        double CSP = col_sum1.dot(col_sum2) / n_col - CF;
        
        double ESP = TSP - GSP - RSP - CSP;
        int DFE = (n_col - 1) * (n_col - 2);
        
        GMP = GSP / (n_gen - 1);
        EMP = ESP / DFE;
        
        if (cov_type == 1) {
          // Genotypic covariance = (GMP - EMP) / t
          covariance = (GMP - EMP) / n_col;
        } else {
          // Phenotypic covariance = Genotypic + EMP
          covariance = (GMP - EMP) / n_col + EMP;
        }
        
      } else {
        // ========== Split Plot Design ==========
        VectorXd main_sum1 = main_sums.col(i);
        VectorXd main_sum2 = main_sums.col(j);
        double MSP = main_sum1.dot(main_sum2) / n_rep - CF;
        
        // Compute main × rep interaction sums
        MatrixXd main_rep_sums(n_main, n_rep);
        main_rep_sums.setZero();
        VectorXd mr_prod1 = VectorXd::Zero(n_main * n_rep);
        VectorXd mr_prod2 = VectorXd::Zero(n_main * n_rep);
        
        for (int obs = 0; obs < n_obs; ++obs) {
          int m = main(obs);
          int r = rep(obs);
          int idx = m * n_rep + r;
          mr_prod1(idx) += trait1(obs);
          mr_prod2(idx) += trait2(obs);
        }
        
        double IMSP = 0.0;
        for (int idx = 0; idx < n_main * n_rep; ++idx) {
          IMSP += mr_prod1(idx) * mr_prod2(idx);
        }
        IMSP = IMSP / (n_obs / (n_main * n_rep)) - CF - MSP - RSP;
        
        double ESP = TSP - GSP - MSP - RSP - IMSP;
        
        int DFE = (n_gen - 1) * (n_main - 1) * (n_rep - 1);
        int DFE_MAIN = (n_main - 1) * (n_rep - 1);
        
        GMP = GSP / (n_gen - 1);
        EMP = ESP / DFE;
        double EMP_MAIN = IMSP / DFE_MAIN;
        
        if (cov_type == 1) {
          // Genotypic covariance = (GMP - EMP) / (r * a)
          covariance = (GMP - EMP) / (n_rep * n_main);
        } else {
          // Phenotypic covariance = Genotypic + EMP
          covariance = (GMP - EMP) / (n_rep * n_main) + EMP;
        }
      }
      
      // Store result (symmetric matrix)
      result(i, j) = covariance;
      if (i != j) {
        result(j, i) = covariance;
      }
    }
  }
  
  return result;
}
