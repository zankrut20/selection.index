# Linear Genomic Eigen Selection Index Methods

## Introduction

This vignette provides an overview of the linear molecular and genomic
eigen selection index methods, which integrate marker trait data and
genome-wide marker values into the calculation. Methods covered include
the Molecular Eigen Selection Index Method (MESIM), Genomic Eigen
Selection Index Method (GESIM), Genome-Wide Linear Eigen Selection Index
Method (GW-ESIM), Restricted Genomic Eigen Selection Index Method
(RGESIM), and Predetermined Proportional Gain Genomic Eigen Selection
Index Method (PPG-GESIM).

## Setup Data Matrices

``` r
library(selection.index)

# Load standard phenotype dataset
data(maize_pheno)

# Define traits and design variables
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
env_col <- "Environment"
genotype_col <- "Genotype"

# Phenotypic variance-covariance matrix (P)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno[[genotype_col]], maize_pheno[[env_col]])

# Genetic variance-covariance matrix (G)
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno[[genotype_col]], maize_pheno[[env_col]])

# For the sake of demonstration within this vignette, we simulate the required molecular/genomic variance components:
set.seed(42)

# Simulate Gamma: Covariance between phenotypes and GEBVs
Gamma <- gmat * 0.85

# Molecular matrices for MESIM
S_M <- gmat * 0.75 # Covariance between phenotypic values and marker scores
S_Mg <- gmat * 0.70 # Covariance between genotypic values and marker scores
S_var <- gmat * 0.80 # Variance-covariance of marker scores

# Genomic matrices for GW-ESIM
G_M <- gmat * 0.82 # Covariance between true genotypic values and marker values
M <- gmat * 0.90 # Variance-covariance matrix of markers
```

## The Molecular Eigen Selection Index Method (MESIM)

The MESIM index integrates marker scores directly:

$$I = \beta\prime_{y}\mathbf{y} + \beta\prime_{s}\mathbf{s} = \begin{bmatrix}
{\beta\prime_{y}} & {\beta\prime_{s}}
\end{bmatrix}\begin{bmatrix}
\mathbf{y} \\
\mathbf{s}
\end{bmatrix} = \beta\prime\mathbf{t}$$

Its optimum coefficients are calculated effectively from the primary
eigenvector resulting from this maximization function:

$$\left( \mathbf{T}^{- 1}\mathbf{\Psi} - \lambda_{M}^{2}\mathbf{I}_{2t} \right)\beta_{M} = \mathbf{0}$$

Using the `mesim` function:

``` r
mes_index <- mesim(pmat, gmat, S_M, S_Mg, S_var)
summary(mes_index)
#> 
#> ==============================================================
#> MOLECULAR EIGEN SELECTION INDEX METHOD (MESIM)
#> Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.1
#> ==============================================================
#> 
#> Selection intensity (k_I): 2.063 
#> Number of traits:          3 
#> 
#> -------------------------------------------------------------
#> INDEX METRICS
#> -------------------------------------------------------------
#>   lambda_M^2 (h^2_I):     1.000798
#>   Accuracy (r_HI):        1.000399
#>   Index Std Dev (sigma_I): 6.315343
#>   Selection Response (R_M): 13.028553
#> 
#> -------------------------------------------------------------
#> PHENOTYPE COEFFICIENTS (b_y)
#> -------------------------------------------------------------
#>           Trait       b_y
#>           Yield  0.000162
#>     PlantHeight -0.010260
#>  DaysToMaturity  0.007350
#> 
#> -------------------------------------------------------------
#> MARKER SCORE COEFFICIENTS (b_s)
#> -------------------------------------------------------------
#>           Trait       b_s
#>           Yield -0.012821
#>     PlantHeight  0.812791
#>  DaysToMaturity -0.582277
#> 
#> -------------------------------------------------------------
#> EXPECTED GENETIC GAINS PER TRAIT (E_M)
#> -------------------------------------------------------------
#>           Trait        E_M
#>           Yield 165.154061
#>     PlantHeight  16.311970
#>  DaysToMaturity  -0.538044
#> 
#> ==============================================================
#> SUMMARY TABLE
#> ==============================================================
#> 
#>     b_y.1    b_y.2   b_y.3     b_s.1    b_s.2     b_s.3      hI2      rHI
#>  0.000162 -0.01026 0.00735 -0.012821 0.812791 -0.582277 1.000798 1.000399
#>   sigma_I      R_M  lambda2
#>  6.315343 13.02855 1.000798
```

## The Linear Genomic Eigen Selection Index Method (GESIM)

The GESIM incorporates Genomic Estimated Breeding Values (GEBVs, denoted
as $\gamma$):

$$I = \beta\prime_{y}\mathbf{y} + \beta\prime_{\gamma}{\mathbf{γ}} = \begin{bmatrix}
{\beta\prime_{y}} & {\beta\prime_{\gamma}}
\end{bmatrix}\begin{bmatrix}
\mathbf{y} \\
{\mathbf{γ}}
\end{bmatrix} = \beta\prime\mathbf{f}$$

The optimum index coefficients are the first eigenvector of:

$$\left( \mathbf{\Phi}^{- 1}\mathbf{A} - \lambda_{G}^{2}\mathbf{I}_{2t} \right)\beta_{G} = \mathbf{0}$$

Using the `gesim` function:

``` r
ges_index <- gesim(pmat, gmat, Gamma)
summary(ges_index)
#> 
#> ==============================================================
#> LINEAR GENOMIC EIGEN SELECTION INDEX METHOD (GESIM)
#> Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.2
#> ==============================================================
#> 
#> Selection intensity (k_I): 2.063 
#> Number of traits:          3 
#> 
#> -------------------------------------------------------------
#> INDEX METRICS
#> -------------------------------------------------------------
#>   lambda_G^2 (h^2_I):     1.000000
#>   Accuracy (r_HI):        1.000000
#>   Index Std Dev (sigma_I): 18221.926505
#>   Selection Response (R_G): 37591.834380
#> 
#> -------------------------------------------------------------
#> PHENOTYPE COEFFICIENTS (b_y)
#> -------------------------------------------------------------
#>           Trait b_y
#>           Yield   0
#>     PlantHeight   0
#>  DaysToMaturity   0
#> 
#> -------------------------------------------------------------
#> GEBV COEFFICIENTS (b_gamma)
#> -------------------------------------------------------------
#>           Trait b_gamma
#>           Yield       1
#>     PlantHeight       0
#>  DaysToMaturity       0
#> 
#> -------------------------------------------------------------
#> EXPECTED GENETIC GAINS PER TRAIT (E_G)
#> -------------------------------------------------------------
#>           Trait        E_G
#>           Yield 37591.8344
#>     PlantHeight   985.3145
#>  DaysToMaturity   547.5454
#> 
#> -------------------------------------------------------------
#> IMPLIED ECONOMIC WEIGHTS (w_G)
#> -------------------------------------------------------------
#>           Trait Implied_w
#>           Yield  0.000255
#>     PlantHeight -0.011029
#>  DaysToMaturity  0.002325
#> 
#> ==============================================================
#> SUMMARY TABLE
#> ==============================================================
#> 
#>  b_y.1 b_y.2 b_y.3 b_gamma.1 b_gamma.2 b_gamma.3 hI2 rHI  sigma_I      R_G
#>      0     0     0         1         0         0   1   1 18221.93 37591.83
#>  lambda2
#>        1
```

## The Genome-Wide Linear Eigen Selection Index Method (GW-ESIM)

When high-density markers covering the whole genome are available, the
GW-ESIM formulation estimates genetic potential uniformly across the
entire marker panel.

Using the `gw_esim` function:

``` r
gw_index <- gw_esim(pmat, gmat, G_M, M)
summary(gw_index)
#> 
#> ==============================================================
#> GENOME-WIDE LINEAR EIGEN SELECTION INDEX METHOD (GW-ESIM)
#> Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.3
#> ==============================================================
#> 
#> Selection intensity (k_I): 2.063 
#> Number of traits:          3 
#> Number of markers:         3 
#> 
#> -------------------------------------------------------------
#> INDEX METRICS
#> -------------------------------------------------------------
#>   lambda_W^2 (h^2_I):     1.000000
#>   Accuracy (r_HI):        1.000000
#>   Index Std Dev (sigma_I): 18750.207685
#>   Selection Response (R_W): 38681.678454
#> 
#> -------------------------------------------------------------
#> PHENOTYPE COEFFICIENTS (b_y)
#> -------------------------------------------------------------
#>           Trait b_y
#>           Yield   0
#>     PlantHeight   0
#>  DaysToMaturity   0
#> 
#> -------------------------------------------------------------
#> MARKER COEFFICIENTS (b_m) - First 10
#> -------------------------------------------------------------
#>  Marker b_m
#>       1   1
#>       2   0
#>       3   0
#> 
#> -------------------------------------------------------------
#> EXPECTED GENETIC GAINS PER TRAIT (E_W)
#> -------------------------------------------------------------
#>           Trait        E_W
#>           Yield 35243.3070
#>     PlantHeight   923.7575
#>  DaysToMaturity   513.3379
#> 
#> ==============================================================
#> SUMMARY TABLE
#> ==============================================================
#> 
#>  n_traits n_markers hI2 rHI  sigma_I      R_W lambda2
#>         3         3   1   1 18750.21 38681.68       1
```

## The Restricted Linear Genomic Eigen Selection Index Method (RGESIM)

RGESIM enables breeders to impose constraints on particular traits such
that their expected genetic advancements are zero, while employing the
genomic capability of the GESIM.

Using the `rgesim` function to restrict improvements to `Ear_Height`
(second trait):

``` r
# Restrict the second trait (PlantHeight)
U_mat <- matrix(c(0, 1, 0), nrow = 1)
rges_index <- rgesim(pmat, gmat, Gamma, U_mat)
summary(rges_index)
#> 
#> ==============================================================
#> RESTRICTED LINEAR GENOMIC EIGEN SELECTION INDEX (RGESIM)
#> Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.4
#> ==============================================================
#> 
#> Selection intensity (k_I): 2.063 
#> Number of traits:          3 
#> Number of restrictions:    1 
#> 
#> -------------------------------------------------------------
#> INDEX METRICS
#> -------------------------------------------------------------
#>   lambda_RG^2 (h^2_I):    1.000000
#>   Accuracy (r_HI):        1.000000
#>   Index Std Dev (sigma_I): 6.999294
#>   Selection Response (R_RG): 14.439543
#> 
#> -------------------------------------------------------------
#> PHENOTYPE COEFFICIENTS (b_y)
#> -------------------------------------------------------------
#>           Trait b_y
#>           Yield   0
#>     PlantHeight   0
#>  DaysToMaturity   0
#> 
#> -------------------------------------------------------------
#> GEBV COEFFICIENTS (b_gamma)
#> -------------------------------------------------------------
#>           Trait   b_gamma
#>           Yield -0.026207
#>     PlantHeight  0.999657
#>  DaysToMaturity  0.000000
#> 
#> -------------------------------------------------------------
#> EXPECTED GENETIC GAINS PER TRAIT (E_RG)
#> -------------------------------------------------------------
#>           Trait       E_RG
#>           Yield -550.97005
#>     PlantHeight    0.00000
#>  DaysToMaturity  -10.82068
#> 
#> -------------------------------------------------------------
#> CONSTRAINT VERIFICATION
#> -------------------------------------------------------------
#> Constrained response (should be near zero):
#> [1] 0
#> 
#> -------------------------------------------------------------
#> IMPLIED ECONOMIC WEIGHTS (w_RG)
#> -------------------------------------------------------------
#>           Trait Implied_w
#>           Yield -0.011032
#>     PlantHeight  0.476755
#>  DaysToMaturity -0.100501
#> 
#> ==============================================================
#> SUMMARY TABLE
#> ==============================================================
#> 
#>  b_y.1 b_y.2 b_y.3 b_gamma.1 b_gamma.2 b_gamma.3 hI2 rHI  sigma_I     R_RG
#>      0     0     0 -0.026207  0.999657         0   1   1 6.999294 14.43954
#>  lambda2
#>        1
```

## The Predetermined Proportional Gain Linear Genomic Eigen Selection Index Method (PPG-GESIM)

This approach aims to acquire customized genetic gains according to a
relative priority proportion assigned distinctly for each trait.

$$\left( \mathbf{T}_{PG} - \lambda_{PG}^{2}\mathbf{I}_{2t} \right)\beta_{PG} = \mathbf{0}$$

Using the `ppg_gesim` function with relative gains heavily weighted
towards `Yield` (trait 3) and `Days_To_Silking` (trait 4):

``` r
# Desired genetic gain proportions: 1 for Yield, 0.5 for DaysToMaturity, 0 for others
d <- c(1, 0, 0.5)
ppg_ges_index <- ppg_gesim(pmat, gmat, Gamma, d)
summary(ppg_ges_index)
#> 
#> ==============================================================
#> PREDETERMINED PROPORTIONAL GAIN GENOMIC EIGEN INDEX (PPG-GESIM)
#> Ceron-Rojas & Crossa (2018) - Chapter 8, Section 8.5
#> ==============================================================
#> 
#> Selection intensity (k_I): 2.063 
#> Number of traits:          3 
#> 
#> -------------------------------------------------------------
#> DESIRED PROPORTIONAL GAINS (d)
#> -------------------------------------------------------------
#>           Trait   d
#>           Yield 1.0
#>     PlantHeight 0.0
#>  DaysToMaturity 0.5
#> 
#> -------------------------------------------------------------
#> INDEX METRICS
#> -------------------------------------------------------------
#>   lambda_PG^2 (h^2_I):    3.500000
#>   Accuracy (r_HI):        1.870829
#>   Index Std Dev (sigma_I): 2056.464205
#>   Selection Response (R_PG): 4242.485656
#> 
#> -------------------------------------------------------------
#> PHENOTYPE COEFFICIENTS (b_y)
#> -------------------------------------------------------------
#>           Trait       b_y
#>           Yield -0.001616
#>     PlantHeight  0.019492
#>  DaysToMaturity -0.610526
#> 
#> -------------------------------------------------------------
#> GEBV COEFFICIENTS (b_gamma)
#> -------------------------------------------------------------
#>           Trait  b_gamma
#>           Yield 0.111810
#>     PlantHeight 0.002665
#>  DaysToMaturity 0.783816
#> 
#> -------------------------------------------------------------
#> EXPECTED GENETIC GAINS PER TRAIT (E_PG)
#> -------------------------------------------------------------
#>           Trait       E_PG
#>           Yield 37151.3513
#>     PlantHeight   973.7697
#>  DaysToMaturity   541.1297
#> 
#> -------------------------------------------------------------
#> PROPORTIONAL GAIN VERIFICATION
#> -------------------------------------------------------------
#> Gain ratios (E_PG / d) - should be approximately constant:
#>           Trait     Ratio
#>           Yield 37151.351
#>     PlantHeight        NA
#>  DaysToMaturity  1082.259
#> Standard deviation of ratios: 25504.7 
#> 
#> -------------------------------------------------------------
#> IMPLIED ECONOMIC WEIGHTS (w_PG)
#> -------------------------------------------------------------
#>           Trait Implied_w
#>           Yield  0.046996
#>     PlantHeight -0.015539
#>  DaysToMaturity  0.004212
#> 
#> ==============================================================
#> SUMMARY TABLE
#> ==============================================================
#> 
#>      b_y.1    b_y.2     b_y.3 b_gamma.1 b_gamma.2 b_gamma.3 hI2      rHI
#>  -0.001616 0.019492 -0.610526   0.11181  0.002665  0.783816 3.5 1.870829
#>   sigma_I     R_PG lambda2
#>  2056.464 4242.486     3.5
```

## Statistical Properties and Efficiency

Just as with phenotypic and standard molecular selection indices, the
efficiency of Genomic Eigen Selection Indices can be evaluated via
accuracy ($\rho_{HI}$), maximized selection response ($R$), and expected
genetic gain per trait ($E$).

For the MESIM index, the parameters are formulated by estimating the
standard deviations of the true genetic merit ($\sigma_{H}$) and the
given index ($\sigma_{I}$). The maximum accuracy obtained through
canonical correlation ($\lambda_{M}$) corresponds directly to the square
root of the primary eigenvalue:

$$\rho_{H_{M}I_{M}} = \frac{\sqrt{\beta\prime_{M}\mathbf{T}_{M}\beta_{M}}}{\sqrt{\beta\prime_{M}\mathbf{T}_{M}\mathbf{\Psi}_{M}^{- 1}\mathbf{T}_{M}\beta_{M}}} = \frac{\sigma_{I_{M}}}{\sigma_{H_{M}}}$$

The selection response is obtained using a standardized selection
intensity factor ($k_{I}$):

$$R_{M} = k_{I}\sqrt{\beta\prime_{M_{1}}\mathbf{T}_{M}\beta_{M_{1}}}$$

Expected genetic gain per trait can then be calculated accordingly:

$$\mathbf{E}_{M} = k_{I}\frac{\mathbf{\Psi}_{M}\beta_{M_{1}}}{\sqrt{\beta\prime_{M_{1}}\mathbf{T}_{M}\beta_{M_{1}}}}$$
