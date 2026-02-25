# Linear Genomic Selection Indices

## Abstract

The linear genomic selection index (LGSI) is a linear combination of
genomic estimated breeding values (GEBVs) used to predict the individual
net genetic merit and select individual candidates from a non-phenotyped
testing population as parents of the next selection cycle. In the LGSI,
phenotypic and marker data from the training population are fitted into
a statistical model to estimate all individual genome marker effects;
these estimates can then be used in subsequent selection cycles to
obtain GEBVs that are predictors of breeding values in a testing
population for which there is only marker information.

Applying the LGSI requires predicting and ranking the net genetic merit
of candidates for selection using only marker information. The combined
LGSI (CLGSI), conversely, uses both phenotypic and GEBV information
jointly to predict the net genetic merit. The CLGSI is typically used in
the training populations to estimate true merit.

This vignette covers calculating LGSI and CLGSI using the
`selection.index` package.

## 5.1 The Linear Genomic Selection Index

### Basic Conditions for Constructing the LGSI

To construct a valid LGSI, four essential conditions must be met: 1. All
marker effects must be estimated simultaneously in the training
population. 2. The estimated marker effects must be used in subsequent
selection cycles to obtain GEBVs that act as predictors of the
individual breeding values in the testing population (where there is
only marker information). 3. The GEBV values must be composed entirely
of the additive genetic effects. 4. Phenotypes must be used to estimate
all marker effects in the training population, and *not* to make
selections in the testing population (Heffner et al., 2009; Lorenz et
al., 2011).

### The LGSI and Its Parameters

The objective of LGSI is to predict the net genetic merit
$H = \mathbf{w}\prime\mathbf{g}$, where $\mathbf{g}$ is a vector of
unobservable true breeding values, and $\mathbf{w}$ is a vector of
economic weights. Let $\mathbf{γ}$ be the vector of known Genomic
Estimated Breeding Values (GEBVs). The basic LGSI equation is:

$$I_{G} = {\mathbf{β}}^{\prime}{\mathbf{γ}}$$

where $\mathbf{β}$ is the unknown variable vector of optimal weights.
For LGSI, it has been shown that maximizing the accuracy of
$\rho_{HI_{G}}$ results simply in ${\mathbf{β}} = \mathbf{w}$.

The main advantage of the LGSI over the conventional Linear Phenotypic
Selection Index (LPSI) lies in the possibility of reducing the intervals
between selection cycles from roughly 4 years to 1.5 years in plant
breeding. The selection response equation adjusts for cycle length
$L_{G}$:

$$R_{I_{G}} = \frac{k_{I}}{L_{G}}\sqrt{\mathbf{w}^{\prime}\mathbf{\Gamma}\mathbf{w}}$$
Where $\mathbf{\Gamma}$ is the covariance matrix of the GEBVs.

### Statistical LGSI Properties

The LGSI operates as a direct genomic extension of the LPSI theory,
retaining the exact same statistical properties: 1. The variance of
$I_{G}$ ($\sigma_{I_{G}}^{2}$) and the covariance between $H$ and
$I_{G}$ ($\sigma_{HI_{G}}$) are equal
($\sigma_{I_{G}}^{2} = \sigma_{HI_{G}}$). 2. The maximized correlation
between $H$ and $I_{G}$ (the LGSI accuracy) is equal to
$\rho_{HI_{G}} = \frac{\sigma_{I_{G}}}{\sigma_{H}}$. 3. The variance of
the predicted error,
$Var\left( H - I_{G} \right) = \left( 1 - \rho_{HI_{G}}^{2} \right)\sigma_{H}^{2}$,
is minimal. 4. The total variance of $H$ explained by $I_{G}$ is
$\sigma_{I_{G}}^{2} = \rho_{HI_{G}}^{2}\sigma_{H}^{2}$.

### LGSI vs LPSI Efficiency

To formalize the efficiency between the Genomic (LGSI) and Phenotypic
(LPSI) selection index efficiency, we use the parameter $\lambda_{G}$:
$$p_{G} = 100\left( \lambda_{G} - 1 \right)$$

If $p_{G} > 0$, LGSI efficiency is greater than LPSI efficiency; if
$p_{G} = 0$, both are equal, and if $p_{G} < 0$, LPSI is more efficient.
Additionally, accounting for the interval lengths between selection
cycles, LGSI will be strictly more efficient than LPSI when:
$$L_{G} < \frac{\rho_{HI_{G}}}{h_{I}}L_{P}$$

### Genomic Estimated Breeding Values (GEBVs)

Because the `selection.index` package does not perform prediction
modeling itself, the GEBVs must be calculated prior. In this example, we
generate GEBVs ($\mathbf{γ}$) using an unpenalized Ridge Regression,
where markers estimate trait values. For a rigorous methodology it is
acceptable to generate GEBVs using Bayesian Alphabets, rrBLUP, or
equivalent prediction models.

``` r
library(selection.index)
library(MASS) # For robust ridge regression implementation (lm.ridge)

# 1. Load the built-in synthetic maize datasets
data("maize_pheno")
data("maize_geno")

# Ensure markers are numeric matrices
X <- as.matrix(maize_geno)

# Ensure pheno lines up perfectly with geno
# (Note: maize_pheno contains 4 repetitions per genotype; we take the means)
Y_agg <- mean_performance(
    data = maize_pheno[, c("Yield", "PlantHeight", "DaysToMaturity")],
    genotypes = maize_pheno$Genotype,
    replications = maize_pheno$Block
)
#> Warning in sqrt(EMS/r): NaNs produced
#> Warning in sqrt(2 * EMS/r): NaNs produced
#> Warning in sqrt(2 * EMS/r): NaNs produced
#> Warning in sqrt(EMS): NaNs produced
#> Warning in sqrt(EMS/r): NaNs produced
#> Warning in sqrt(2 * EMS/r): NaNs produced
#> Warning in sqrt(2 * EMS/r): NaNs produced
#> Warning in sqrt(EMS): NaNs produced
#> Warning in sqrt(EMS/r): NaNs produced
#> Warning in sqrt(2 * EMS/r): NaNs produced
#> Warning in sqrt(2 * EMS/r): NaNs produced
#> Warning in sqrt(EMS): NaNs produced

# Sort both datasets to ensure identical ordering
Y_agg <- Y_agg[order(Y_agg$Genotypes), ]
X <- X[order(rownames(X)), ]

# Ensure overlapping genotypes
common_genotypes <- intersect(Y_agg$Genotypes, rownames(X))
Y_agg <- Y_agg[Y_agg$Genotypes %in% common_genotypes, ]
X <- X[rownames(X) %in% common_genotypes, ]

# Select only multi-traits relevant to breeding
# Yield, PlantHeight, DaysToMaturity
Y <- as.matrix(Y_agg[, c("Yield", "PlantHeight", "DaysToMaturity")])
```

We now calculate $\mathbf{γ}$ (the `gebv_mat`):

``` r
# Simulate Genomic Estimated Breeding Values (GEBVs) using Ridge Regression
# In best practice, you'd use cross-validation or a separate training/testing set
# We use lambda = 100 to handle the p >> n dimensionality problem
gebv_mat <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
colnames(gebv_mat) <- colnames(Y)

for (i in seq_len(ncol(Y))) {
    # Fit a ridge regression model for trait `i` using markers `X`
    model_ridge <- lm.ridge(Y[, i] ~ X, lambda = 100)

    # Predict values using coef() which correctly un-scales the coefficients: Include intercept
    betas <- coef(model_ridge)
    intercept <- betas[1]
    beta <- betas[-1] # Exclude intercept

    # Calculate the GEBV for the trait
    gebv_mat[, i] <- intercept + (X %*% beta)
}

head(gebv_mat, 3)
#>         Yield PlantHeight DaysToMaturity
#> [1,] 8241.151    219.4007       115.8738
#> [2,] 7965.285    208.8754       116.4623
#> [3,] 7832.597    204.5494       116.4592
```

### Implementing LGSI

With `gebv_mat` in hand, it is simple to calculate the LGSI metrics
using
[`lgsi()`](https://zankrut20.github.io/selection.index/reference/lgsi.md):

``` r
# 2. Compute Trait Covariances
pmat <- cov(Y) # Phenotypic Covariance (Approximation)

# In optimal practice, the true Genomic Covariance Matrix (\Gamma) is
# estimated using Restricted Maximum Likelihood (REML) utilizing both general
# phenotypic and genotypic information. Here, we simulate \Gamma as proportional
# to the phenotypic variance (assuming a high heritability correlation).
gmat <- pmat * 0.4 # Approximate Genotypic Covariance (Approximating \Gamma)

# 3. Define Economic Weights
weights <- data.frame(
    Trait = c("Yield", "PlantHeight", "DaysToMaturity"),
    Weight = c(5, -0.1, -0.1)
)
wmat <- weight_mat(weights)

# 4. Calculate Linear Genomic Selection Index (LGSI)
# For the testing population where we only use genomic values
lgsi_result <- lgsi(
    gebv_mat = gebv_mat,
    gmat = gmat,
    wmat = wmat
)

# Output Summary
lgsi_result$summary
#>   b.1  b.2  b.3       GA      PRE hI2 rHI
#> 1   5 -0.1 -0.1 1718.838 171883.8   1   1
```

The LGSI function output shows the calculated standard deviation
($\sigma_{I}$), index heritability ($h^{2}$), the overall net genetic
advance ($GA$), and expected gains per trait ($\Delta H$) associated
with genomic selection.

------------------------------------------------------------------------

## 5.2 The Combined Linear Genomic Selection Index

The Combined LGSI (CLGSI) developed by Dekkers (2007) is a slightly
modified version of the LMSI, which, instead of using the raw marker
scores, uses the GEBVs and the phenotypic information jointly to predict
the net genetic merit.

The CLGSI model represents:
$$I_{C} = {\mathbf{β}}_{\mathbf{y}}^{\prime}\mathbf{y} + {\mathbf{β}}_{G}^{\prime}{\mathbf{γ}}$$

Where $\mathbf{y}$ is the predicted phenotype and $\mathbf{γ}$ is the
subset of genetic markers. As $\left. p\rightarrow\infty \right.$, the
true variance approaches the genomic predicted variance, and CLGSI
becomes identical to LGSI. Because the CLGSI uses GEBVs and phenotypic
information jointly, it can ideally complement prediction workflows
inside training populations where field trials are simultaneously held.

### Statistical Properties of the CLGSI

Similar to LGSI, the CLGSI holds strict conditions ensuring theoretical
validation constraints: 1. $\sigma_{I_{C}}^{2} = \sigma_{HI_{C}}$,
meaning the variance of $I_{C}$ and the covariance between $H$ and
$I_{C}$ are unified. 2. The maximized correlation between $H$ and
$I_{C}$ operates strictly as
$\rho_{HI_{C}} = \frac{\sigma_{I_{C}}}{\sigma_{H}}$. 3. The variance of
the predicted error,
$Var\left( H - I_{C} \right) = \left( 1 - \rho_{HI_{C}}^{2} \right)\sigma_{H}^{2}$,
is mathematically minimal. 4. The total variance of $H$ explained by
$I_{C}$ scales linearly via
$\sigma_{I_{C}}^{2} = \rho_{HI_{C}}^{2}\sigma_{H}^{2}$.

### Implementing CLGSI

You can measure the CLGSI utilizing the
[`clgsi()`](https://zankrut20.github.io/selection.index/reference/clgsi.md)
function and feeding it both the observed phenotypes and estimated
genomic variants.

``` r
clgsi_result <- clgsi(
    phen_mat = Y, # Observed phenotypic data
    gebv_mat = gebv_mat, # Genomic Estimated Breeding Values
    pmat = pmat, # Expected Phenotypic traits covariance
    gmat = gmat, # Expected Genotypic traits covariance
    wmat = wmat
)

clgsi_result$summary
#>       b_y.1    b_y.2    b_y.3    b_g.1    b_g.2    b_g.3       GA      PRE hI2
#> 1 -241.2056 -69.0278 -11711.7 301.3417 -12.6428 13421.05 9325.477 932547.7   1
#>   rHI
#> 1   1
```

Notice that the Combined LGSI has two sets of optimized unitless $\beta$
weightings—one for the phenotype values (`b_y`), and one for the overall
marker value (`b_g`) per trait. Selection Response ($R$) represents the
final expected genetic aggregate improvement.

## Limitations and Caveats

1.  **Dimensionality Problem**: Often $p \gg n$, requiring
    regularization techniques like Ridge Regression (as demonstrated in
    `lm.ridge`) to derive the Genomic Estimated Breeding Values without
    singular/non-invertible covariance anomalies.
2.  **Cycle Length**: The genomic methods (LGSI, CLGSI) may occasionally
    show a lower index accuracy unadjusted than the phenotypic methods
    initially, but typically provide vastly more improvement per unit of
    time ($L_{G}$) because of how heavily compressed cycle iterations
    are when phenotyping is skipped.
