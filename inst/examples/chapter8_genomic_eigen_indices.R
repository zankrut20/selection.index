# ==============================================================================
# Chapter 8: Linear Molecular and Genomic Eigen Selection Index Methods
# ==============================================================================
#
# This script demonstrates all five genomic eigen selection index methods
# from Chapter 8 of Cerón-Rojas & Crossa (2018).
#
# Methods covered:
# 1. MESIM    - Molecular Eigen Selection Index Method
# 2. GESIM    - Linear Genomic Eigen Selection Index Method
# 3. GW-ESIM  - Genome-Wide Linear Eigen Selection Index Method
# 4. RGESIM   - Restricted Linear Genomic Eigen Selection Index Method
# 5. PPG-GESIM- Predetermined Proportional Gain Genomic Eigen Selection Index
#
# These methods extend the eigen-based approach (Chapter 7) to genomic/molecular
# data, maximizing accuracy without requiring pre-specified economic weights.
#
# ==============================================================================

library(selection.index)

# Load example data
data(seldata)

# Compute variance-covariance matrices
gmat <- gen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])
pmat <- phen_varcov(seldata[, 3:9], seldata[, 2], seldata[, 1])

# Number of traits
n_traits <- nrow(gmat)
trait_names <- colnames(gmat)

cat("\n=============================================================================\n")
cat("Chapter 8: Genomic Eigen Selection Index Methods\n")
cat("=============================================================================\n\n")

cat("Number of traits:", n_traits, "\n")
cat("Trait names:", paste(trait_names, collapse = ", "), "\n\n")

# ==============================================================================
# 8.1  MESIM – Molecular Eigen Selection Index Method
# ==============================================================================

cat("\n--- 8.1 MESIM: Molecular Eigen Selection Index Method ---\n\n")

# Simulate marker score covariance matrices
# In practice, these would be computed from actual marker scores and genetic data
set.seed(123)

# S_M: Cov(y, s) - Covariance between phenotypes and marker scores
# Used in T_M (phenotypic variance matrix)
S_M <- gmat * 0.7
S_M <- (S_M + t(S_M)) / 2  # Ensure symmetry

# S_Mg: Cov(g, s) - Covariance between genetic values and marker scores  
# Used in Ψ_M (genetic variance matrix)
# Theoretically distinct from S_M, but often S_Mg ≈ S_M when Cov(e,s) ≈ 0
S_Mg <- gmat * 0.65  # Slightly lower to reflect pure genetic correlation
S_Mg <- (S_Mg + t(S_Mg)) / 2  # Ensure symmetry

# S_var: Var(s) - Variance of marker scores
# Used in both T_M and Ψ_M (assuming no measurement error in scores)
S_var <- gmat * 0.8
S_var <- (S_var + t(S_var)) / 2  # Ensure symmetry

# Fit MESIM with all parameters (most theoretically rigorous)
result_mesim <- mesim(pmat, gmat, S_M, S_Mg = S_Mg, S_var = S_var, 
                       selection_intensity = 2.063)

# Alternative usages:
# Standard (assumes Cov(e,s) ≈ 0, so Cov(g,s) ≈ Cov(y,s)):
#   result_mesim <- mesim(pmat, gmat, S_M, S_var = S_var)
# Backward compatible (Chapter 8.1 simplified notation):
#   result_mesim <- mesim(pmat, gmat, S_M)

cat("MESIM Summary:\n")
print(result_mesim$summary)

cat("\nPhenotype coefficients (b_y):\n")
print(round(result_mesim$b_y, 4))

cat("\nMarker score coefficients (b_s):\n")
print(round(result_mesim$b_s, 4))

cat("\nExpected genetic gains per trait (E_M):\n")
print(round(result_mesim$E_M, 4))

cat("\nIndex heritability (h²_I):", round(result_mesim$hI2, 4), "\n")
cat("Index accuracy (r_HI):", round(result_mesim$rHI, 4), "\n")
cat("Selection response (R_M):", round(result_mesim$R_M, 4), "\n")

# ==============================================================================
# 8.2  GESIM – Linear Genomic Eigen Selection Index Method
# ==============================================================================

cat("\n\n--- 8.2 GESIM: Linear Genomic Eigen Selection Index Method ---\n\n")

# Simulate GEBV covariance matrix
# In practice, this would be computed from GEBVs obtained via genomic prediction
# Gamma represents Cov(phenotypes, GEBVs)
set.seed(456)
Gamma <- gmat * 0.8  # Assume 80% GEBV-phenotype covariance
Gamma <- (Gamma + t(Gamma)) / 2  # Ensure symmetry

# Fit GESIM
result_gesim <- gesim(pmat, gmat, Gamma, selection_intensity = 2.063)

cat("GESIM Summary:\n")
print(result_gesim$summary)

cat("\nPhenotype coefficients (b_y):\n")
print(round(result_gesim$b_y, 4))

cat("\nGEBV coefficients (b_gamma):\n")
print(round(result_gesim$b_gamma, 4))

cat("\nExpected genetic gains per trait (E_G):\n")
print(round(result_gesim$E_G, 4))

cat("\nImplied economic weights:\n")
print(round(result_gesim$implied_w, 4))

cat("\nIndex heritability (h²_I):", round(result_gesim$hI2, 4), "\n")
cat("Index accuracy (r_HI):", round(result_gesim$rHI, 4), "\n")
cat("Selection response (R_G):", round(result_gesim$R_G, 4), "\n")

# ==============================================================================
# 8.3  GW-ESIM – Genome-Wide Linear Eigen Selection Index Method
# ==============================================================================

cat("\n\n--- 8.3 GW-ESIM: Genome-Wide Linear Eigen Selection Index Method ---\n\n")

# Simulate genome-wide marker data
# In practice, these would be computed from actual marker scores
set.seed(789)
N_markers <- 50  # Use modest number for demonstration

# G_M: Covariance between phenotypes and marker scores (n_traits x N_markers)
G_M <- matrix(rnorm(n_traits * N_markers, mean = 0, sd = 0.5), 
              nrow = n_traits, ncol = N_markers)

# M: Variance-covariance matrix of marker scores (N_markers x N_markers)
M <- diag(N_markers) + matrix(rnorm(N_markers^2, sd = 0.1), N_markers, N_markers)
M <- (M + t(M)) / 2  # Make symmetric
M <- M + diag(0.1, N_markers)  # Ensure positive definite

# Fit GW-ESIM
result_gw_esim <- gw_esim(pmat, gmat, G_M, M, selection_intensity = 2.063)

cat("GW-ESIM Summary:\n")
print(result_gw_esim$summary)

cat("\nPhenotype coefficients (b_y):\n")
print(round(result_gw_esim$b_y, 4))

cat("\nFirst 10 marker coefficients (b_m):\n")
print(round(head(result_gw_esim$b_m, 10), 4))

cat("\nExpected genetic gains per trait (E_W):\n")
print(round(result_gw_esim$E_W, 4))

cat("\nIndex heritability (h²_I):", round(result_gw_esim$hI2, 4), "\n")
cat("Index accuracy (r_HI):", round(result_gw_esim$rHI, 4), "\n")
cat("Selection response (R_W):", round(result_gw_esim$R_W, 4), "\n")

# ==============================================================================
# 8.4  RGESIM – Restricted Linear Genomic Eigen Selection Index Method
# ==============================================================================

cat("\n\n--- 8.4 RGESIM: Restricted Linear Genomic Eigen Selection Index Method ---\n\n")

# Define restrictions: hold first trait at zero gain
# This demonstrates restricting gains on specific traits while optimizing others
# IMPORTANT: U_mat must be (n_traits x n_restrictions)
# The 2r restriction rule is automatically applied by rgesim()
U_mat <- matrix(0, nrow = n_traits, ncol = 1)
U_mat[1, 1] <- 1  # Restrict trait 1

cat("Restriction matrix U (holding Trait 1 at zero gain):\n")
print(U_mat)

# Fit RGESIM
result_rgesim <- rgesim(pmat, gmat, Gamma, U_mat, selection_intensity = 2.063)

cat("\nRGESIM Summary:\n")
print(result_rgesim$summary)

cat("\nPhenotype coefficients (b_y):\n")
print(round(result_rgesim$b_y, 4))

cat("\nGEBV coefficients (b_gamma):\n")
print(round(result_rgesim$b_gamma, 4))

cat("\nExpected genetic gains per trait (E_RG):\n")
print(round(result_rgesim$E_RG, 4))

cat("\nConstrained response (should be near zero):\n")
print(round(result_rgesim$constrained_response, 6))

cat("\nImplied economic weights:\n")
print(round(result_rgesim$implied_w, 4))

cat("\nIndex heritability (h²_I):", round(result_rgesim$hI2, 4), "\n")
cat("Index accuracy (r_HI):", round(result_rgesim$rHI, 4), "\n")
cat("Selection response (R_RG):", round(result_rgesim$R_RG, 4), "\n")

# Verify constraint satisfaction
cat("\nVerification: Trait 1 gain should be near zero:\n")
cat("E_RG[1] =", round(result_rgesim$E_RG[1], 6), "\n")

# ==============================================================================
# 8.5  PPG-GESIM – Predetermined Proportional Gain Genomic Eigen Selection Index
# ==============================================================================

cat("\n\n--- 8.5 PPG-GESIM: Predetermined Proportional Gain Genomic Eigen Selection Index ---\n\n")

# Define desired proportional gains
# Example: Want gains in ratio 2:1:3 for first three traits, equal for rest
d <- c(2, 1, 3, 1, 1, 1, 1)

cat("Desired proportional gains (d):\n")
print(d)

# Fit PPG-GESIM
result_ppg_gesim <- ppg_gesim(pmat, gmat, Gamma, d, selection_intensity = 2.063)

cat("\nPPG-GESIM Summary:\n")
print(result_ppg_gesim$summary)

cat("\nPhenotype coefficients (b_y):\n")
print(round(result_ppg_gesim$b_y, 4))

cat("\nGEBV coefficients (b_gamma):\n")
print(round(result_ppg_gesim$b_gamma, 4))

cat("\nExpected genetic gains per trait (E_PG):\n")
print(round(result_ppg_gesim$E_PG, 4))

cat("\nGain ratios (E_PG / d) - should be approximately constant:\n")
print(round(result_ppg_gesim$gain_ratios, 4))

cat("\nImplied economic weights:\n")
print(round(result_ppg_gesim$implied_w, 4))

cat("\nIndex heritability (h²_I):", round(result_ppg_gesim$hI2, 4), "\n")
cat("Index accuracy (r_HI):", round(result_ppg_gesim$rHI, 4), "\n")
cat("Selection response (R_PG):", round(result_ppg_gesim$R_PG, 4), "\n")

# Verify proportional gains
cat("\nVerification: Check if gain ratios are constant:\n")
gain_ratios_sd <- sd(result_ppg_gesim$gain_ratios, na.rm = TRUE)
cat("Standard deviation of gain ratios:", round(gain_ratios_sd, 6), "\n")
cat("(Small value indicates proportional gains are achieved)\n")

# ==============================================================================
# Comparison of Methods
# ==============================================================================

cat("\n\n=============================================================================\n")
cat("COMPARISON OF ALL CHAPTER 8 METHODS\n")
cat("=============================================================================\n\n")

comparison_df <- data.frame(
  Method = c("MESIM", "GESIM", "GW-ESIM", "RGESIM", "PPG-GESIM"),
  h2_I = c(result_mesim$hI2, result_gesim$hI2, result_gw_esim$hI2, 
           result_rgesim$hI2, result_ppg_gesim$hI2),
  r_HI = c(result_mesim$rHI, result_gesim$rHI, result_gw_esim$rHI,
           result_rgesim$rHI, result_ppg_gesim$rHI),
  R = c(result_mesim$R_M, result_gesim$R_G, result_gw_esim$R_W,
        result_rgesim$R_RG, result_ppg_gesim$R_PG),
  lambda2 = c(result_mesim$lambda2, result_gesim$lambda2, result_gw_esim$lambda2,
              result_rgesim$lambda2, result_ppg_gesim$lambda2),
  Constraints = c("None", "None", "None", "Zero gain on Trait 1", 
                  "Proportional gains")
)

cat("Summary comparison:\n")
print(round(comparison_df[, -6], 4))

cat("\n\nMethod characteristics:\n")
cat("  MESIM    : Uses marker scores directly\n")
cat("  GESIM    : Uses GEBVs, provides implied weights\n")
cat("  GW-ESIM  : Genome-wide approach with all markers\n")
cat("  RGESIM   : Restricts specific traits to zero gain\n")
cat("  PPG-GESIM: Enforces predetermined proportional gains\n")

cat("\n=============================================================================\n")
cat("END OF CHAPTER 8 EXAMPLES\n")
cat("=============================================================================\n")
