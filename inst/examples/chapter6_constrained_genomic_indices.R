# ==============================================================================
# EXAMPLES: Constrained Linear Genomic Selection Indices (Chapter 6)
# ==============================================================================
# 
# This script demonstrates the four constrained genomic selection indices
# implemented from Chapter 6 of Cerón-Rojas & Crossa (2018).
#
# All functions follow the package conventions:
# - Use C++ primitives from math_primitives.cpp for computations
# - All selection logic implemented in R
# - Modular design in constrained_genomic_indices.R
#
# ==============================================================================

library(selection.index)

# Load example data
data(seldata)

# Compute variance-covariance matrices
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Display trait names
cat("Traits:", colnames(gmat), "\n\n")

# Economic weights (example from weight dataset)
data(weight)
wmat <- weight_mat(weight)
w <- wmat[, 1]  # Use first column

# ==============================================================================
# EXAMPLE 1: RLGSI (Restricted Linear Genomic Selection Index)
# ==============================================================================
# Scenario: Restrict traits 2 (EHT) and 4 (ED) to zero gain
# 
# Mathematical form (Eq. 6.3):
#   β_RG = K_G w  where U'Γβ_RG = 0
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("EXAMPLE 1: RLGSI - Restrict traits 2 and 4 to zero gain\n")
cat(strrep("=", 78), "\n\n")

# Use gmat as Gamma (GEBV variance-covariance)
Gamma <- gmat

result_rlgsi <- rlgsi(
  Gamma = Gamma,
  wmat = w,
  restricted_traits = c(2, 4),  # Restrict EHT and ED
  k_I = 2.063,  # 10% selection intensity
  L_G = 1
)

cat("RLGSI Coefficients:\n")
print(round(result_rlgsi$b, 4))

cat("\nExpected Genetic Gains (E_RG):\n")
print(round(result_rlgsi$E, 4))

cat("\nSelection Response (R_RG):", round(result_rlgsi$R, 4), "\n")

cat("\nConstraint Check (should be ~0 for restricted traits):\n")
print(round(result_rlgsi$E[c(2, 4)], 6))

cat("\nSummary:\n")
print(result_rlgsi$summary)

cat("\n")

# ==============================================================================
# EXAMPLE 2: PPG-LGSI (Predetermined Proportional Gains)
# ==============================================================================
# Scenario: Achieve proportional gains with ratio 3:2:1:0.5:0:1:2
# 
# Mathematical form (Eq. 6.9):
#   β_PG = β_RG + θ_G U(U'ΓU)^(-1)d
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("EXAMPLE 2: PPG-LGSI - Predetermined proportional gains\n")
cat(strrep("=", 78), "\n\n")

# Desired proportional gains
d <- c(3, 2, 1, 0.5, 0, 1, 2)

result_ppg_lgsi <- ppg_lgsi(
  Gamma = Gamma,
  d = d,
  wmat = w,
  k_I = 2.063,
  L_G = 1
)

cat("PPG-LGSI Coefficients:\n")
print(round(result_ppg_lgsi$b, 4))

cat("\nExpected Genetic Gains (E_PG):\n")
print(round(result_ppg_lgsi$E, 4))

cat("\nDesired Proportional Gains (d):\n")
print(d)

cat("\nGain Ratios (E / d) - should be constant:\n")
print(round(result_ppg_lgsi$gain_ratios, 4))

cat("\nProportionality constant (θ_G):", round(result_ppg_lgsi$theta_G, 4), "\n")

cat("\nSummary:\n")
print(result_ppg_lgsi$summary)

cat("\n")

# ==============================================================================
# EXAMPLE 3: CRLGSI (Combined Restricted Linear Genomic Selection Index)
# ==============================================================================
# Scenario: Combine phenotypes + GEBVs with restrictions
# 
# Mathematical form (Eq. 6.13):
#   β_CR = K_C β_C  where constraints on genetic gains
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("EXAMPLE 3: CRLGSI - Combined data with restrictions\n")
cat(strrep("=", 78), "\n\n")

# Simulate phenotypes and GEBVs for demonstration
set.seed(123)
n_genotypes <- 100
n_traits <- ncol(gmat)

# Simulate phenotypes (with environmental noise)
phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                   nrow = n_genotypes, ncol = n_traits)
colnames(phen_mat) <- colnames(gmat)

# Simulate GEBVs (correlated with genetic component)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
                   nrow = n_genotypes, ncol = n_traits)
colnames(gebv_mat) <- colnames(gmat)

result_crlgsi <- crlgsi(
  phen_mat = phen_mat,
  gebv_mat = gebv_mat,
  pmat = pmat,
  gmat = gmat,
  wmat = w,
  restricted_traits = c(2, 4),  # Restrict EHT and ED
  reliability = 0.7,  # 70% reliability of GEBVs
  k_I = 2.063,
  L_I = 1
)

cat("CRLGSI Coefficients for Phenotypes (b_y):\n")
print(round(result_crlgsi$b_y, 4))

cat("\nCRLGSI Coefficients for GEBVs (b_g):\n")
print(round(result_crlgsi$b_g, 4))

cat("\nExpected Genetic Gains (E_CR):\n")
print(round(result_crlgsi$E, 4))

cat("\nSelection Response (R_CR):", round(result_crlgsi$R, 4), "\n")

cat("\nSummary:\n")
print(result_crlgsi$summary)

cat("\n")

# ==============================================================================
# EXAMPLE 4: CPPG-LGSI (Combined Predetermined Proportional Gains)
# ==============================================================================
# Scenario: Most general - combine phenotypes + GEBVs with proportional gains
# 
# Mathematical form (Eq. 6.16):
#   β_CP = β_CR + θ_CP δ_CP
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("EXAMPLE 4: CPPG-LGSI - Combined data with proportional gains\n")
cat(strrep("=", 78), "\n\n")

result_cppg_lgsi <- cppg_lgsi(
  phen_mat = phen_mat,
  gebv_mat = gebv_mat,
  pmat = pmat,
  gmat = gmat,
  d = d,  # Same proportional gains as Example 2
  wmat = w,
  reliability = 0.7,
  k_I = 2.063,
  L_I = 1
)

cat("CPPG-LGSI Coefficients for Phenotypes (b_y):\n")
print(round(result_cppg_lgsi$b_y, 4))

cat("\nCPPG-LGSI Coefficients for GEBVs (b_g):\n")
print(round(result_cppg_lgsi$b_g, 4))

cat("\nExpected Genetic Gains (E_CP):\n")
print(round(result_cppg_lgsi$E, 4))

cat("\nDesired Proportional Gains (d):\n")
print(d)

cat("\nGain Ratios (E / d):\n")
print(round(result_cppg_lgsi$gain_ratios, 4))

cat("\nProportionality constant (θ_CP):", round(result_cppg_lgsi$theta_CP, 4), "\n")

cat("\nSummary:\n")
print(result_cppg_lgsi$summary)

cat("\n")

# ==============================================================================
# COMPARISON SUMMARY
# ==============================================================================

cat(strrep("=", 78), "\n")
cat("COMPARISON OF METHODS\n")
cat(strrep("=", 78), "\n\n")

comparison <- data.frame(
  Method = c("RLGSI", "PPG-LGSI", "CRLGSI", "CPPG-LGSI"),
  Data_Type = c("GEBV only", "GEBV only", "Phenotype + GEBV", "Phenotype + GEBV"),
  Constraint = c("Zero gain", "Proportional", "Zero gain", "Proportional"),
  Response = round(c(result_rlgsi$R, result_ppg_lgsi$R, 
                     result_crlgsi$R, result_cppg_lgsi$R), 4),
  Accuracy = round(c(result_rlgsi$rHI, result_ppg_lgsi$rHI,
                     result_crlgsi$rHI, result_cppg_lgsi$rHI), 4)
)

print(comparison)

cat("\n")
cat("Notes:\n")
cat("- Combined methods (CRLGSI, CPPG-LGSI) may show higher accuracy\n")
cat("- PPG methods achieve proportional gains as specified\n")
cat("- Restricted methods maintain zero gain for constrained traits\n")
cat("- All calculations use C++ primitives from math_primitives.cpp\n")
cat("- All logic implemented in R for transparency\n")

# ==============================================================================
# END OF EXAMPLES
# ==============================================================================

cat("\n")
cat(strrep("=", 78), "\n")
cat("Examples completed successfully!\n")
cat(strrep("=", 78), "\n")
