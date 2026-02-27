# Generate synthetic maize data for the selection.index package

set.seed(2026)

# 1. Phenotypic Data (maize_pheno)
# Simulating 100 genotypes across 2 environments and 3 blocks for 3 traits
n_genotypes <- 100
n_env <- 2
n_blocks <- 3

maize_pheno <- expand.grid(
  Genotype = paste0("G", sprintf("%03d", 1:n_genotypes)),
  Environment = paste0("Env", 1:n_env),
  Block = paste0("B", 1:n_blocks)
)

# Simulate traits with some correlation
base_yield <- rnorm(nrow(maize_pheno), mean = 8000, sd = 500)
base_height <- base_yield * 0.02 + rnorm(nrow(maize_pheno), mean = 50, sd = 10)
base_maturity <- 120 - (base_yield * 0.001) + rnorm(nrow(maize_pheno), mean = 5, sd = 2)

maize_pheno$Yield <- round(base_yield, 2)
maize_pheno$PlantHeight <- round(base_height, 2)
maize_pheno$DaysToMaturity <- round(base_maturity, 0)

# 2. Genomic Data (maize_geno)
# Simulating a matrix of 100 genotypes x 500 SNP markers (coded 0, 1, 2)
n_markers <- 500
maize_geno <- matrix(
  sample(0:2, n_genotypes * n_markers, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
  nrow = n_genotypes,
  ncol = n_markers
)
rownames(maize_geno) <- paste0("G", sprintf("%03d", 1:n_genotypes))
colnames(maize_geno) <- paste0("SNP_", sprintf("%03d", 1:n_markers))

# Save the generated data to the data/ directory
usethis::use_data(maize_pheno, maize_geno, overwrite = TRUE)
