#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

//' C++ ANOVA Statistics Iterator
//'
//' @description
//' Optimized C++ implementation for computing ANOVA statistics across multiple traits.
//' Replaces R loop in mean_performance() that calls design_stats() for each trait.
//' Processes all traits in a single pass using vectorized operations.
//'
//' @param data_mat Numeric matrix of trait observations (n_obs x n_traits)
//' @param gen_idx Integer vector of genotype indices (1-based converted to 0-based internally)
//' @param rep_idx Integer vector of replication indices (1-based converted to 0-based internally)
//' @param col_idx Integer vector of column indices for LSD (NULL for RCBD/SPD)
//' @param main_idx Integer vector of main plot indices for SPD (NULL for RCBD/LSD)
//' @param design_type Integer: 1=RCBD, 2=LSD, 3=SPD
//'
//' @return List containing matrices:
//'   - GMS: Genotype Mean Squares (1 x n_traits)
//'   - EMS: Error Mean Squares (1 x n_traits)
//'   - EMS_MAIN: Main plot error mean squares for SPD (1 x n_traits)
//'   - DFG: Degrees of freedom for genotypes (scalar)
//'   - DFE: Degrees of freedom for error (scalar)
//'   - DFE_MAIN: Degrees of freedom for main plot error in SPD (scalar)
//'   - n_rep: Number of replications (scalar)
//'   - n_gen: Number of genotypes (scalar)
//'   - n_main: Number of main plots for SPD (scalar)
//'
//' @details
//' Algorithm optimizations:
//' - Single pass through data for all traits (eliminates R loop overhead)
//' - Pre-computes all grouped sums for reuse across traits
//' - Vectorized operations using Eigen
//' - Returns statistics needed for F-tests, CV%, heritability, etc.
//'
//' @keywords internal
// [[Rcpp::export]]
List cpp_anova_iterator(
    const Eigen::Map<Eigen::MatrixXd>& data_mat,
    const Eigen::Map<Eigen::VectorXi>& gen_idx,
    const Eigen::Map<Eigen::VectorXi>& rep_idx,
    Nullable<Eigen::Map<Eigen::VectorXi>> col_idx = R_NilValue,
    Nullable<Eigen::Map<Eigen::VectorXi>> main_idx = R_NilValue,
    int design_type = 1  // 1=RCBD, 2=LSD, 3=SPD
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
  
  // Pre-allocate result vectors for all traits
  VectorXd GMS_vec(n_traits);
  VectorXd EMS_vec(n_traits);
  VectorXd EMS_MAIN_vec(n_traits);
  
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
  
  // Process each trait
  for (int t = 0; t < n_traits; ++t) {
    VectorXd trait = data_mat.col(t);
    
    // Correction Factor
    double CF = (grand_totals(t) * grand_totals(t)) / n_obs;
    
    // Total Sum of Squares
    double TSS = trait.dot(trait) - CF;
    
    // Genotype Sum of Squares
    VectorXd gen_sum = gen_sums.col(t);
    double GSS = gen_sum.dot(gen_sum) / n_rep - CF;
    
    // Replication Sum of Squares
    VectorXd rep_sum = rep_sums.col(t);
    double RSS = rep_sum.dot(rep_sum) / n_gen - CF;
    
    int DFG = n_gen - 1;
    int DFR = n_rep - 1;
    int DFE, DFE_MAIN = 0;
    double ESS, ESS_MAIN = 0;
    
    if (design_type == 1) {
      // ========== RCBD ==========
      ESS = TSS - GSS - RSS;
      DFE = DFG * DFR;
      
      GMS_vec(t) = GSS / DFG;
      EMS_vec(t) = ESS / DFE;
      EMS_MAIN_vec(t) = 0.0;  // Not applicable for RCBD
      
    } else if (design_type == 2) {
      // ========== Latin Square Design ==========
      VectorXd col_sum = col_sums.col(t);
      double CSS = col_sum.dot(col_sum) / n_col - CF;
      
      ESS = TSS - GSS - RSS - CSS;
      DFE = (n_col - 1) * (n_col - 2);
      
      GMS_vec(t) = GSS / DFG;
      EMS_vec(t) = ESS / DFE;
      EMS_MAIN_vec(t) = 0.0;  // Not applicable for LSD
      
    } else {
      // ========== Split Plot Design ==========
      VectorXd main_sum = main_sums.col(t);
      double MSS = main_sum.dot(main_sum) / n_rep - CF;
      
      // Compute main × rep interaction sums
      VectorXd mr_prod = VectorXd::Zero(n_main * n_rep);
      
      for (int obs = 0; obs < n_obs; ++obs) {
        int m = main(obs);
        int r = rep(obs);
        int idx = m * n_rep + r;
        mr_prod(idx) += trait(obs);
      }
      
      double IMSS = 0.0;
      for (int idx = 0; idx < n_main * n_rep; ++idx) {
        IMSS += mr_prod(idx) * mr_prod(idx);
      }
      IMSS = IMSS / (n_obs / (n_main * n_rep)) - CF - MSS - RSS;
      
      // Sub-plot (genotype) sum of squares
      ESS = TSS - GSS - MSS - RSS - IMSS;
      
      // Degrees of freedom
      int DFM = n_main - 1;
      DFE = (n_gen - 1) * (n_main - 1) * (n_rep - 1);
      DFE_MAIN = (n_main - 1) * (n_rep - 1);
      ESS_MAIN = IMSS;
      
      GMS_vec(t) = GSS / DFG;
      EMS_vec(t) = ESS / DFE;
      EMS_MAIN_vec(t) = ESS_MAIN / DFE_MAIN;
    }
  }
  
  // Calculate degrees of freedom (same for all traits)
  int DFG = n_gen - 1;
  int DFE, DFE_MAIN = 0;
  
  if (design_type == 1) {
    DFE = DFG * (n_rep - 1);
  } else if (design_type == 2) {
    DFE = (n_col - 1) * (n_col - 2);
  } else {
    DFE = (n_gen - 1) * (n_main - 1) * (n_rep - 1);
    DFE_MAIN = (n_main - 1) * (n_rep - 1);
  }
  
  return List::create(
    Named("GMS") = GMS_vec,
    Named("EMS") = EMS_vec,
    Named("EMS_MAIN") = EMS_MAIN_vec,
    Named("DFG") = DFG,
    Named("DFE") = DFE,
    Named("DFE_MAIN") = DFE_MAIN,
    Named("n_rep") = n_rep,
    Named("n_gen") = n_gen,
    Named("n_main") = n_main
  );
}
