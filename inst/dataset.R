set.seed(42)

# ==============================================================================
# 1. PARAMETERS
# ==============================================================================
n_genotypes <- 150
n_reps <- 3
n_traits <- 4
trait_names <- c("Yield", "PlantHeight", "EarHeight", "FloweringTime")

# True Genetic Correlation Matrix (High correlation between heights)
Sigma_G <- matrix(c(
  20,  10,   8,  -5,
  10,  40,  25,   2,
  8,  25,  30,   1,
  -5,   2,   1,  15
), nrow=4, ncol=4, dimnames = list(trait_names, trait_names))

# Environmental Variance (Heritability ~ 0.4 to 0.6)
Sigma_E <- Sigma_G * 1.5

# ==============================================================================
# 2. GENERATE TRUE BREEDING VALUES (TBV)
# ==============================================================================
# Multivariate normal simulation
TBV <- MASS::mvrnorm(n = n_genotypes, mu = rep(100, n_traits), Sigma = Sigma_G)
rownames(TBV) <- paste0("G", 1:n_genotypes)

# ==============================================================================
# 3. CREATE PHENOTYPIC DATASET (maize_pheno) - For Chapters 2, 3
# ==============================================================================
# RCBD Design: Genotypes x Reps
pheno_data <- data.frame()

for (r in 1:n_reps) {
  # Add environmental error specific to this rep
  error <- MASS::mvrnorm(n = n_genotypes, mu = rep(0, n_traits), Sigma = Sigma_E)

  # Phenotype = TBV + Error
  Y <- TBV + error

  tmp <- data.frame(
    Genotype = rownames(TBV),
    Rep = r,
    Y
  )
  pheno_data <- rbind(pheno_data, tmp)
}

maize_pheno <- pheno_data
head(maize_pheno)

# ==============================================================================
# 4. CREATE MARKER DATASET (maize_markers) - For Chapter 4
# ==============================================================================
n_markers <- 500
# Simulate markers (-1, 0, 1 coding)
markers <- matrix(sample(c(-1, 0, 1), n_genotypes * n_markers, replace = TRUE),
                  nrow = n_genotypes, ncol = n_markers)
rownames(markers) <- rownames(TBV)
colnames(markers) <- paste0("M", 1:n_markers)

maize_markers <- markers

# ==============================================================================
# 5. CREATE GENOMIC ESTIMATED BREEDING VALUES (maize_genomics) - For Ch 5,6,8,9
# ==============================================================================
# Simulate GEBVs by adding some noise to True Breeding Values (Accuracy ~ 0.7)
# GEBV = TBV + PredictionError
Sigma_GEBV_Error <- Sigma_G * 0.4
GEBV_error <- MASS::mvrnorm(n = n_genotypes, mu = rep(0, n_traits), Sigma = Sigma_GEBV_Error)

maize_genomics <- TBV + GEBV_error
rownames(maize_genomics) <- rownames(TBV)

# ==============================================================================
# 6. ECONOMIC WEIGHTS (maize_weights)
# ==============================================================================
maize_weights <- c(10, 2, 1, -2) # High Yield, moderate Height, negative Flowering
names(maize_weights) <- trait_names

library(usethis)
use_data(maize_pheno, maize_genomics, maize_markers, maize_weights, overwrite = TRUE)
