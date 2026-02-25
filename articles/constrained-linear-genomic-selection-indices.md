# Constrained Linear Genomic Selection Indices

## Introduction

This vignette demonstrates the methodology behind **Constrained Linear
Genomic Selection Indices (LGSI)** based on Chapter 6 of the book
*Linear Selection Indices in Modern Plant Breeding*. The objective of
constrained selection indices is to maximize genetic advance for
specific traits while imposing defined restrictions on other traits.

In a genomic selection context, these indices use Genomic Estimated
Breeding Values (GEBVs). The core indices covered here include:

1.  **Restricted Linear Genomic Selection Index (RLGSI)**: Aims to
    improve specific traits while leaving others fixed (null
    restriction).
2.  **Predetermined Proportional Gain Linear Genomic Selection Index
    (PPG-LGSI)**: Aims to achieve specific, predetermined proportions of
    change in the population mean of designated traits.
3.  **Combined RLGSI (CRLGSI)**: Extends RLGSI to use phenotypic and
    GEBV information jointly.
4.  **Combined PPG-LGSI (CPPG-LGSI)**: Extends PPG-LGSI to use
    phenotypic and GEBV information jointly.

## Environment Setup and Data Preparation

We will utilize synthetic phenotypic and genotypic maize data provided
in the `selection.index` package to illustrate these methods.

### 1. Estimating GEBVs

Before applying these indices, we first estimate genomic breeding
values. We use Ridge Regression (`lm.ridge`) from the `MASS` package as
a rapid and simple method to simulate this step linking marker data to
phenotype data.

``` r
library(selection.index)
library(MASS) # For Ridge regression

# Load the maize datasets
data("maize_pheno", package = "selection.index")
data("maize_geno", package = "selection.index")

# Extract only the traits we care about to prevent missing value logic failures
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
maize_pheno_clean <- na.omit(maize_pheno[, c("Genotype", "Block", traits)])

# Calculate mean performance for phenotypic data
# Providing only the traits to `data` to accurately generate the table
pheno_means_full <- mean_performance(
    data = maize_pheno_clean[, traits],
    genotypes = maize_pheno_clean$Genotype,
    replications = maize_pheno_clean$Block,
    method = "Mean"
)

# Remove the trailing 9 rows composed of summary text-stats automatically generated
pheno_means <- pheno_means_full[1:(nrow(pheno_means_full) - 9), ]

# Extract numerical phenotype matrix
Y <- as.matrix(pheno_means[, traits])
mode(Y) <- "numeric"
rownames(Y) <- pheno_means$Genotypes

# Format genotype matrix
X <- as.matrix(maize_geno)

# Ensure matching dimensions and align the matrices
common_genotypes <- intersect(rownames(Y), rownames(X))

# Filter and sort matrices to match genotypes sequence exactly
Y_filtered <- Y[common_genotypes, ]
X_filtered <- X[common_genotypes, ]

# Calculate GEBVs using a simple Ridge Regression model per trait
gebvs_sim <- matrix(0, nrow = nrow(Y_filtered), ncol = ncol(Y_filtered))
colnames(gebvs_sim) <- colnames(Y_filtered)
rownames(gebvs_sim) <- rownames(Y_filtered)
lambda_ridge <- 0.1

for (j in 1:ncol(Y_filtered)) {
    # Fit ridge regression to simulate marker effects
    model_ridge <- lm.ridge(Y_filtered[, j] ~ X_filtered, lambda = lambda_ridge)
    beta_ridge <- coef(model_ridge)[-1]

    # Predict GEBVs
    gebvs_sim[, j] <- X_filtered %*% matrix(beta_ridge, ncol = 1)
}

# Define Covariance Matrices
Gamma <- cov(gebvs_sim) # Genomic Covariance Matrix estimated from simulated GEBVs
P_matrix <- cov(Y_filtered) # Phenotypic Covariance Matrix
C_matrix <- Gamma # Genetic Covariance (Assuming Gamma captures genetic additive cov)

# Define vector of economic weights for Yield, PlantHeight, DaysToMaturity
w <- matrix(c(5, -0.1, -0.1), ncol = 1)
```

## 6.1 The Restricted Linear Genomic Selection Index (RLGSI)

The RLGSI is applied in a testing population when only GEBVs are
available. The goal is to maximize the response for some traits while
restricting the expected generic advance of other traits to exactly
zero.

### Mathematical Definition

The expected genetic advance for RLGSI is restricted by:
$Cov\left( I_{RG},\mathbf{U}\prime\mathbf{g} \right) = \mathbf{U}\prime\mathbf{\Gamma}{\mathbf{β}}_{RG} = \mathbf{0}$,
where $\mathbf{U}\prime$ is a matrix of 1s and 0s identifying the
restricted traits.

Following the maximization step of the selection response, the index
coefficients ${\mathbf{β}}_{RG}$ optimally resolve as:
$${\mathbf{β}}_{RG} = \mathbf{K}_{G}\mathbf{w}$$ where
$\mathbf{K}_{G} = \left\lbrack \mathbf{I} - \mathbf{Q}_{G} \right\rbrack$
and
$\mathbf{Q}_{G} = \mathbf{U}(\mathbf{U}\prime\mathbf{\Gamma}\mathbf{U})^{- 1}\mathbf{U}\prime\mathbf{\Gamma}$.

### RLGSI Implementation Example

Suppose we want to maximize the gain using our three traits, but we want
to hold **Plant Height** and **Days To Maturity** completely fixed ($0$
gain).

``` r
# Define the restriction matrix U
# We restrict traits 2 (PlantHeight) and 3 (DaysToMaturity).
# Trait 1 (Yield) is left unrestricted.
U_rlgsi <- matrix(c(
    0, 0, # Yield unrestricted
    1, 0, # PlantHeight restricted
    0, 1 # DaysToMaturity restricted
), nrow = 3, byrow = TRUE)

# Calculate RLGSI
rlgsi_res <- rlgsi(Gamma = Gamma, wmat = w, U = U_rlgsi)

# View the resulting index coefficients and Expected Genetic Gain
cat("RLGSI Coefficients (beta_RG):\n")
#> RLGSI Coefficients (beta_RG):
print(rlgsi_res$b)
#> [1]    5.0000 -111.1516  175.9708

cat("\nExpected Genetic Gain per Trait:\n")
#> 
#> Expected Genetic Gain per Trait:
print(rlgsi_res$Summary[, "Delta_H", drop = FALSE])
#> NULL
```

Notice in the output that the Expected Genetic Gain (`Delta_H`) for
traits 2 and 3 evaluates to a precision limit of $0$, fulfilling the
null-gain restriction.

## 6.2 The Predetermined Proportional Gain Linear Genomic Selection Index (PPG-LGSI)

Instead of forcing a trait to remain fixed at $0$, we might want a
trait’s mean value to undergo a highly specific, predefined amount of
selection shift.

### Mathematical Definition

Let
$\mathbf{d}\prime = \left\lbrack d_{1},d_{2},\ldots,d_{r} \right\rbrack$
be the proportional gain values requested by the breeder. Similar to the
RLGSI, the coefficient array evaluates identically with a targeted
constant proportionality factor restricting expected proportional
shifts: $${\mathbf{β}}_{PG} = \mathbf{K}_{P}\mathbf{w}$$

### PPG-LGSI Implementation Example

Suppose we want to target **Yield** to proportionally increase by $7.0$
units, and identically target **Plant Height** to slightly decrease by
relatively $- 3.0$ proportion targets.

``` r
# Define the restriction matrix U for PPG
# Restrict traits 1 (Yield) and 2 (PlantHeight) to predetermined scalars
U_ppg <- matrix(c(
    1, 0, # Yield restricted
    0, 1, # PlantHeight restricted
    0, 0 # DaysToMaturity unrestricted
), nrow = 3, byrow = TRUE)

# Define the predetermined values vector 'd'
d_vec <- matrix(c(7.0, -3.0), ncol = 1)

# Calculate PPG-LGSI
ppg_res <- ppg_lgsi(Gamma = Gamma, d = d_vec, wmat = w, U = U_ppg)

cat("PPG-LGSI Coefficients (beta_PG):\n")
#> PPG-LGSI Coefficients (beta_PG):
print(ppg_res$b)
#> [1]   0.2441 -11.1971  -0.1000
```

## 6.3 Combined Restricted Linear Genomic Selection Index (CRLGSI)

The CRLGSI extends the basic RLGSI by integrating *both* phenotypic
information and genomic breeding values structurally as blocks in a
combined index $I_{C}$. It is mainly applied to a training population
where joint data is uniformly available.

Its mathematical structure models parallel coefficient distributions
across phenotypic and genomic parameters uniformly:
$${\mathbf{β}}_{CR} = \mathbf{K}_{C}{\mathbf{β}}_{C}$$

For an index combining phenotypes and GEBVs, restrictions must
concurrently be assigned to both empirical dimensions globally.

``` r
# 1. Provide Combined Matrices
# The combined Phenotypic-Genomic covariance matrix T_C
T_C <- rbind(
    cbind(P_matrix, Gamma),
    cbind(Gamma, Gamma)
)

# The combined Genetic-Genomic covariance matrix Psi_C
Psi_C <- rbind(
    cbind(C_matrix, Gamma),
    cbind(Gamma, Gamma)
)

# 2. Define Restriction Matrix U_C
# Note: For combined indices of `t` traits, U_C spans exactly 2t rows.
# Restricting Trait 1 (Yield) on both empirical and genomic parameters:
U_crlgsi <- matrix(0, nrow = 6, ncol = 2)
U_crlgsi[1, 1] <- 1 # Trait 1 phenotypic component
U_crlgsi[4, 2] <- 1 # Trait 1 genomic component

# 3. Calculate CRLGSI
crlgsi_res <- crlgsi(T_C = T_C, Psi_C = Psi_C, wmat = w, U = U_crlgsi)

cat("CRLGSI Combined Coefficients (beta_CR):\n")
#> CRLGSI Combined Coefficients (beta_CR):
print(crlgsi_res$b)
#> [1] -6.843307e-05 -1.079212e-03 -5.103276e-02  2.017002e-03 -9.892085e-02
#> [6] -4.895374e-02
```

## 6.4 Combined Predetermined Proportional Gain Linear Genomic Selection Index (CPPG-LGSI)

The CPPG-LGSI combines the dynamic proportional gain vectors to the
$T_{C}$ and $\Psi_{C}$ data blocks identically leveraging the
[`cppg_lgsi()`](https://zankrut20.github.io/selection.index/reference/cppg_lgsi.md)
logic.

### CPPG-LGSI Implementation Example

Suppose we arbitrarily evaluate **Yield** to predetermined restriction
constants.

``` r
# Target Yield using predetermined combined proportions d
d_combined <- matrix(c(7.0, 3.5), ncol = 1)

# Calculate CPPG-LGSI dynamically
cppg_res <- cppg_lgsi(T_C = T_C, Psi_C = Psi_C, d = d_combined, wmat = w, U = U_crlgsi)

cat("CPPG-LGSI Coefficients (beta_CP):\n")
#> CPPG-LGSI Coefficients (beta_CP):
print(cppg_res$b)
#> [1]  4.491869e-06  7.083960e-05  3.349690e-03  4.998047e+00 -7.083568e-05
#> [6] -3.350577e-03
```

## Summary Comparison

Below we extract the standard metrics associated with each index model
(i.e Selection Response $R$):

``` r
comparison_df <- data.frame(
    Index = c("RLGSI", "PPG-LGSI", "CRLGSI", "CPPG-LGSI"),
    Selection_Response = c(
        rlgsi_res$R,
        ppg_res$R,
        crlgsi_res$R,
        cppg_res$R
    ),
    Overall_Genetic_Advance = c(
        rlgsi_res$GA,
        ppg_res$GA,
        crlgsi_res$GA,
        cppg_res$GA
    )
)

print(comparison_df)
#>       Index Selection_Response Overall_Genetic_Advance
#> 1     RLGSI       1490.8413294            1490.8413294
#> 2  PPG-LGSI        102.6128081             102.9592598
#> 3    CRLGSI          0.9390596               0.9390601
#> 4 CPPG-LGSI       2123.2848118            2123.2848118
```
