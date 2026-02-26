# Linear Phenotypic Eigen Selection Index Methods

## Introduction

Based on canonical correlation, singular value decomposition (SVD), and
linear phenotypic selection indices theory, the **Eigen Selection Index
Method (ESIM)**, the **Restricted ESIM (RESIM)**, and the
**Predetermined Proportional Gain ESIM (PPG-ESIM)** use only phenotypic
information to predict the net genetic merit and select superior
genotypes.

While the Linear Phenotypic Selection Index (LPSI) assumes economic
weights are known, the ESIM family estimates them in each selection
cycle to maximize the correlation between the index and the net genetic
merit (i.e., accuracy) and the selection response. ESIM is unrestricted,
whereas RESIM fixes specific traits to zero genetic gain, and PPG-ESIM
ensures traits gain in strict prefixed proportions.

This vignette details the theoretical foundation of these three eigen
analysis indices alongside practical workflows for estimating them using
the `selection.index` package.

### Sample Data Preparation

Throughout this vignette, we will reuse the `maize_pheno` package
dataset to calculate our phenotypic ($\mathbf{P}$) and genotypic
($\mathbf{G}$ or $\mathbf{C}$) covariance matrices.

``` r
library(selection.index)

# Load the built-in maize phenotype dataset
data("maize_pheno")

# Extract the traits of interest
traits <- c("Yield", "PlantHeight", "DaysToMaturity")

# Calculate Genotypic (gmat) and Phenotypic (pmat) covariance matrices
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)
```

## 1. The Linear Phenotypic Eigen Selection Index Method (ESIM)

The theoretical selection response of the ESIM can be written as:

$$R_{I} = k_{I}\sigma_{H}\rho_{HI},\quad(7.1)$$

where $k_{I}$ is the selection intensity, $\sigma_{H}$ is the standard
deviation of the genetic merit, and $\rho_{HI}$ is the correlation
between the trait index and the unobservable genetic merit.

To maximize the index selection response without pre-defined economic
weights, ESIM translates the index estimation into a generalized
eigenvalue problem:

$$\left( \mathbf{P}^{- 1}\mathbf{C} - \lambda_{E}^{2}\mathbf{I} \right)\mathbf{b}_{E} = \mathbf{0},\quad(7.10)$$

where $\mathbf{P}^{- 1}\mathbf{C}$ is the multi-trait heritability
matrix. The optimal ESIM vector of coefficients $\mathbf{b}_{E}$ is the
first eigenvector of $\mathbf{P}^{- 1}\mathbf{C}$, and its accuracy is
the square root of the first eigenvalue ($\lambda_{E}$).

The maximized ESIM selection response can then be written as:

$$R_{E} = k_{I}\sqrt{\mathbf{b}\prime_{E}\mathbf{P}\mathbf{b}_{E}}\quad(7.15)$$

### Practical Application: ESIM

In the package, `esim` maximizes the genetic response without
pre-defined weights for `pmat` and `gmat`.

``` r
# Compute the linear phenotypic eigen selection index
esim_res <- esim(
  pmat = pmat,
  gmat = gmat,
  selection_intensity = 2.063
)

# View the summary
print(esim_res$summary)
#>   Index  lambda2      hI2      rHI  sigma_I Delta_G      b.1      b.2      b.3
#> 1     1 1771.794 1771.794 42.09268 3.286805 6.78068 0.002336 0.009329 0.999954

# View expected genetic gains per trait (Delta_G)
print(esim_res$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>      700636.74       18363.40       10206.29
```

## 2. The Restricted Eigen Selection Index Method (RESIM)

Similar to the RLPSI, the RESIM fixes certain traits so that their
expected genetic gain remains null (zero). Instead of modifying known
economic weights under restrictions, RESIM computes the restricted
eigenvalues natively.

When an orthogonal set of restrictions $\mathbf{U}\prime$ is imposed,
the optimized eigenvalue equation becomes:

$$\left( \mathbf{K}\mathbf{P}^{- 1}\mathbf{C} - \lambda_{R}^{2}\mathbf{I}_{t} \right)\mathbf{b}_{R} = \mathbf{0}.\quad(7.31)$$

where $\mathbf{K} = \lbrack\mathbf{I} - \mathbf{Q}\rbrack$,
$\mathbf{Q} = \mathbf{P}^{- 1}\mathbf{C}\mathbf{U}\left( \mathbf{U}\prime\mathbf{C}\mathbf{P}^{- 1}\mathbf{C}\mathbf{U} \right)^{- 1}\mathbf{U}\prime\mathbf{C}$.
The vector $\mathbf{b}_{R}$ is the restriction-corrected eigenvector.

### Practical Application: RESIM

We use `resim` and restrict `PlantHeight` (trait \#2) to have strictly
zero expected genetic gain.

``` r
# We restrict PlantHeight (the second trait in our data frame)
resim_res <- resim(
  pmat = pmat,
  gmat = gmat,
  restricted_traits = c(2),
  selection_intensity = 2.063
)

# Expected genetic gains per trait
print(resim_res$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>  -1.057493e+02   1.324103e-07  -2.276986e+00
```

Notice that the expected genetic gain for our restricted trait is
mathematically forced to zero.

## 3. The Predetermined Proportional Gain ESIM (PPG-ESIM)

Instead of null gains, the PPG-ESIM allows the breeder to pre-set
optimal levels of growth proportion on specific traits while maintaining
the eigen maximization properties on the non-restricted traits.

This changes the eigen mapping matrix using the proportional
restrictions Mallard operator $\mathbf{M}\prime$:

$$\left( \mathbf{K}_{P}\mathbf{P}^{- 1}\mathbf{C} - \lambda_{P}^{2}\mathbf{I}_{t} \right)\mathbf{b}_{P} = \mathbf{0},\quad(7.43)$$

where $\mathbf{K}_{P}$ is derived from the predetermined proportional
vector parameter $\mathbf{d}$.

### Practical Application: PPG-ESIM

Let’s assume we want a strict $5\%$ gain in PlantHeight and a $10\%$
reduction in DaysToMaturity proportionately compared to Yield.

``` r
# Provide the vector d of desired predetermined proportional gains
# The indices represent non-zero targets. For exactly corresponding target magnitudes, we provide d.
d_vector <- c(0, 0.5, -1) # e.g. proportional mapping

# We use the ppg_esim function
ppgesim_res <- ppg_esim(
  pmat = pmat,
  gmat = gmat,
  d = d_vector,
  selection_intensity = 2.063
)

# View the proportionality of the genetic gains
print(ppgesim_res$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>    -2597.62024      -68.25086      -37.50434
```

Through this mechanism, breeders manipulate genetic covariance paths
analytically accurately.

## 4. Statistical Properties of Eigen Selection Indices

### Variance of Predicted Error (VPE)

The precision of the index in predicting the net genetic merit is
evaluated by its Variance of Predicted Error:

$$VPE\left( I_{E} \right) = \sigma_{H_{E}}^{2} - \sigma_{I_{E}}^{2}\quad(7.19)$$

The relative mean squared effect of $I_{E}$ in predicting $H_{E}$ is
equivalent to its reliability:

$$\frac{Var\left( H_{E} \right) - VPE\left( I_{E} \right)}{Var\left( H_{E} \right)} = \rho_{H_{E}I_{E}}^{2}\quad(7.20)$$

where $\rho_{H_{E}I_{E}}^{2}$ is the squared accuracy of the eigen
index.

### Canonical Correlation Connection

The ESIM optimization maps directly to Canonical Correlation Theory
(Hotelling 1936). The objective is to find a linear combination of
traits having maximal correlation with the net genetic merit. In the
ESIM context, $\mathbf{b}_{E}$ is exactly the first canonical vector of
the matrix $\mathbf{P}^{- 1}\mathbf{C}$, and its associated accuracy
represents the highest achievable canonical correlation between the
phenotypic indices and the genotypic values.
