# Linear Marker and Genome-Wide Selection Indices

## Introduction

There are two main linear selection indices employed in marker-assisted
selection (MAS) to predict the net genetic merit and to select
individual candidates as parents for the next generation: the **linear
marker selection index (LMSI)** and the **genome-wide LMSI (GW-LMSI)**.

Both indices maximize the expected genetic gain per trait by combining
molecular and phenotypic data. The LMSI requires a prior step of
selecting markers significantly linked to the quantitative trait loci
(QTL) to calculate a “marker score”. The GW-LMSI is a single-stage
procedure that treats information at each individual marker as a
separate trait, entering all markers and phenotypic information into a
single index.

In this vignette, we demonstrate how to apply these indices using the
`selection.index` package on phenotypic (`maize_pheno`) and genotypic
(`maize_geno`) data.

### Sample Data Preparation

We use the synthetic maize datasets built into the package to calculate
the phenotypic, genotypic, marker, and covariance matrices necessary for
molecular index modeling.

``` r
library(selection.index)

# Load the built-in maize phenotype and genotype datasets
data("maize_pheno")
data("maize_geno")

# 1. Prepare Phenotypic Data
# We select three traits for our demonstration
traits <- c("Yield", "PlantHeight", "DaysToMaturity")
phen_mat <- as.matrix(maize_pheno[, traits])

# Calculate Genotypic (gmat) and Phenotypic (pmat) covariance matrices
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)

# Aggregate phenotypic data by Genotype to match marker data dimensions
agg_pheno <- aggregate(maize_pheno[, traits], by = list(Genotype = maize_pheno$Genotype), FUN = mean)
phen_mat <- as.matrix(agg_pheno[, traits])

# 2. Prepare Genotypic Data
# Extract only the marker columns (dropping the Genotype ID column)
marker_mat <- as.matrix(maize_geno[, -1])

# 3. Define Economic Weights
# Suppose we want to improve Yield strongly, decrease Plant Height, and lightly increase Days to Maturity
wmat <- weight_mat(data.frame(Trait = traits, Weight = c(2, -1, 0.5)))
```

## 1. The Linear Marker Selection Index (LMSI)

The Linear Marker Selection Index (LMSI) requires an explicit prior
step: calculating the **marker score** ($\mathbf{s}$).

Let $\mathbf{y}_{i} = \mathbf{g}_{i} + \mathbf{e}_{i}$ be the vector of
phenotypic traits. The unobservable genetic value $g_{i}$ can be
predicted using a linear combination of the $M$ selected markers
structurally linked to the QTL that affect the $i$th trait:

$$\mathbf{s}_{i} = \sum\limits_{j = 1}^{M}{\mathbf{θ}}_{j}\mathbf{x}_{j}$$

Where $\mathbf{s}_{i}$ is the marker score, ${\mathbf{θ}}_{j}$ is the
regression coefficient, and $\mathbf{x}_{j}$ is the coded value of the
$j$th markers.

The LMSI combines the phenotypic scores ($\mathbf{y}$) and marker scores
($\mathbf{s}$) to predict the net genetic merit $H$:

\$\$
I_M={\boldsymbol{\beta}}\_y^{\prime}\mathbf{y}+{\boldsymbol{\beta}}\_s^{\prime}\mathbf{s}=\left\[{\boldsymbol{\beta}}\_y^{\prime}\kern0.5em
{\boldsymbol{\beta}}\_s^{\prime}\right\]\left\[\begin{array}{c}\mathbf{y}\\
{}\mathbf{s}\end{array}\right\] \$\$

### Practical Application: LMSI

To apply the LMSI, we must first estimate the marker scores for each
trait. A simple way is to perform multiple linear regression of the
phenotypic values on the marker data, selecting only significant
markers, and using those coefficients to calculate $\mathbf{s}$.

For simplicity in this tutorial, we will use all markers in a linear
model (ridge regression) to predict the traits.

``` r
# Calculate marker scores via a simple Ridge Regression
# (We use MASS::lm.ridge or glmnet in practice, but for the vignette we approximate
# the prediction using the marker matrix directly)
library(MASS)

marker_scores <- matrix(0, nrow = nrow(phen_mat), ncol = ncol(phen_mat))
colnames(marker_scores) <- colnames(phen_mat)

for (i in seq_len(ncol(phen_mat))) {
  # Fit ridge regression to handle multicollinearity among markers
  # Note: A lambda must be chosen carefully in real scenarios
  fit <- lm.ridge(phen_mat[, i] ~ marker_mat, lambda = 10)
  # Calculate predicted marker scores
  marker_scores[, i] <- scale(marker_mat, center = fit$xm, scale = fit$scales) %*% fit$coef + fit$ym
}
```

Now that we have $\mathbf{p}\mathbf{m}\mathbf{a}\mathbf{t}$,
$\mathbf{g}\mathbf{m}\mathbf{a}\mathbf{t}$, the phenotypic matrix
$\mathbf{y}$, and the estimated marker score matrix $\mathbf{s}$, we can
calculate the **Linear Marker Selection Index**:

``` r
# Calculate the LMSI
lmsi_res <- lmsi(
  phen_mat = phen_mat,
  marker_scores = marker_scores,
  pmat = pmat,
  gmat = gmat,
  wmat = wmat
)

# View the LMSI coefficient summary
print(lmsi_res$summary)
#>            Trait   Component            b    w   Delta_H
#> 1          Yield   Phenotype     521.4950  2.0 731752.99
#> 2    PlantHeight   Phenotype    2085.7769 -1.0  19178.95
#> 3 DaysToMaturity   Phenotype  223294.4615  0.5  10659.56
#> 4          Yield MarkerScore    -532.8818  2.0 731752.99
#> 5    PlantHeight MarkerScore   -2003.1664 -1.0  19178.95
#> 6 DaysToMaturity MarkerScore -229077.9808  0.5  10659.56

# Display the expected genetic gains per trait
print(lmsi_res$Delta_G)
#> NULL
```

The output yields the relative contribution weights of the phenotypes
(`b_y.*`) versus the markers (`b_s.*`), and calculates the overall
Expected Genetic Advance.

## 2. The Genome-Wide Linear Marker Selection Index (GW-LMSI)

The genome-wide linear selection index (GW-LMSI) bypasses the two-stage
regression scoring step of the LMSI. Instead, it treats information at
each individual marker as a separate trait. All $m$ marker information
enters together with phenotypic information into a single index.

The GW-LMSI combines the trait phenotypic values ($\mathbf{y}$) and the
molecular information ($\mathbf{m}$) to predict $H$:

\$\$
I_W={\boldsymbol{\beta}}\_y^{\prime}\mathbf{y}+{\boldsymbol{\beta}}\_m^{\prime}\mathbf{m}=\left\[{\boldsymbol{\beta}}\_y^{\prime}\kern0.5em
{\boldsymbol{\beta}}\_m^{\prime}\right\]\left\[\begin{array}{c}\mathbf{y}\\
{}\mathbf{m}\end{array}\right\] \$\$

Because the number of markers $m$ often wildly exceeds the number of
individuals $n$ (the $p \gg n$ problem), the covariance matrix of the
markers ($\mathbf{M}$) and the covariance between phenotypes and markers
($\mathbf{W}$) are generally singular.

Thus, calculating the index coefficients ${\mathbf{β}}_{y}$ and
${\mathbf{β}}_{m}$ typically requires generalized inverses or
regularization (e.g. Ridge Regression / Tikhonov regularization applied
to the covariance matrices).

### Practical Application: GW-LMSI

The `selection.index` package manages the complex construction of the
$\mathbf{M}$ and $\mathbf{W}$ covariance matrices asymptotically and
handles singularity issues using Tikhonov regularization (`lambda`).

``` r
# We use all 500 markers in the maize_geno matrix
# Since p (500) approaches n (600), we supply a ridge regularization lambda
# to ensure the covariance matrices are invertible.
gw_lmsi_res <- gw_lmsi(
  marker_mat = marker_mat,
  trait_mat = phen_mat,
  gmat = gmat,
  wmat = wmat,
  lambda = 0.05
)

# Display the phenotypic weighting coefficients
print(gw_lmsi_res$b_y)
#> NULL

# We can also check the condition number to observe the matrix stability
print(gw_lmsi_res$condition_number)
#> [1] 6.629131

# And summarize the overall index statistics
print(gw_lmsi_res$summary)
#>      Metric    Value
#> 1       rHI   0.0252
#> 2       hI2   1.0000
#> 3   sigma_I 397.3752
#> 4         R 819.7850
#> 5        GA 840.7368
#> 6       PRE       NA
#> 7 n_markers 499.0000
#> 8    lambda   0.0500
```

The GW-LMSI automatically models all the molecular variance in one step.
While it requires computationally heavier matrix inversion (often
mitigated by regularization), it tends to predict the net genetic merit
more efficiently than standard LMSI when trait heritability is very low
or when capturing many small QTL regressions requires immense sample
sizes.

## 3. Advanced Concepts and Limitations

### 3.1 Basic Conditions for Constructing the LMSI

According to Lande and Thompson (1990), applying the LMSI correctly
requires several primary statistical and biological conditions to be met
in the breeding population:

1.  **Linkage Disequilibrium:** The markers and the quantitative trait
    loci (QTL) must be in linkage disequilibrium (LD).
2.  **Additivity:** QTL effects should be combined additively both
    within and between loci (no epistasis).
3.  **Coupling Mode:** The QTL should be in coupling mode where parental
    lines hold uniformly positive or negative alleles.
4.  **Genetic Architecture:** Traits should be affected by a few QTLs
    with large effects rather than exclusively by many small-effect
    QTLs.
5.  **Heritability:** The traits must inherently have a low heritability
    ($h^{2}$). If heritabilities are high, phenotypic selection alone is
    immensely effective.
6.  **Marker Identification:** Markers correlated with the traits of
    interest must first be successfully identified.

### 3.2 The “LMSI Paradox” and Efficiency Comparisons

Will a linear marker selection index always dramatically outperform
standard phenotypic selection? Not necessarily.

Lande and Thompson (1990) showed that LMSI efficiency relative to basic
phenotypic selection ($\lambda_{M}$) for a single trait is:

$$\lambda_{M} = \frac{R_{M}}{R} = \sqrt{\frac{q}{h^{2}} + \frac{(1 - q)^{2}}{1 - {qh}^{2}}}$$

Where $q = \frac{\sigma_{s}^{2}}{\sigma_{g}^{2}}$ is the proportion of
additive genetic variance definitively explained by the markers, and
$h^{2}$ is the trait heritability.

This equation demonstrates the **“LMSI Paradox”** mathematically defined
by Knapp (1998): as $h^{2}$ trends toward zero, $\lambda_{M}$ approaches
infinity. Thus, the massive relative advantage of marker-assisted
selection over phenotypic selection occurs explicitly when phenotypic
trait heritability is extremely low, and the markers successfully
capture a massive fraction of genetic variance. Otherwise, the LMSI
offers diminishing returns and can fix QTLs at an overly fast rate.

### 3.3 Statistical Properties of Molecular Indices

Assuming $H$ and the index $I_{M}$ share a bivariate joint normal
distribution, the statistical properties mirror the standard phenotypic
LPSI:

1.  The variance of the index equals the covariance between the net
    genetic merit and the index: $\sigma_{I_{M}}^{2} = \sigma_{HI_{M}}$.
2.  The maximized correlation (Accuracy) is
    $\rho_{HI_{M}} = \frac{\sigma_{I_{M}}}{\sigma_{H}}$.
3.  The Variance of the Predicted Error
    $Var\left( H - I_{M} \right) = \left( 1 - \rho_{HI_{M}}^{2} \right)\sigma_{H}^{2}$
    is uniquely minimized.

Because molecular indices incorporate *both* the phenotypic and
genotypic marker data, the LMSI accuracy will technically always be
higher than the standard LPSI accuracy.

### 3.4 Advanced Matrix Estimation Methods

While the
[`gw_lmsi()`](https://zankrut20.github.io/selection.index/reference/gw_lmsi.md)
and
[`lmsi()`](https://zankrut20.github.io/selection.index/reference/lmsi.md)
functions in the `selection.index` package abstract away the underlying
algebra, molecular index computation routinely suffers from the
$p \gg n$ covariance problem (vastly more markers than recorded
individuals).

Estimating the marker score variances ($\mathbf{S}$) and the massive
$\mathbf{M}$ and $\mathbf{W}$ covariance matrices using least-squares
often leads to: 1. Marker covariance values erroneously overestimating
the true genetic covariance ($\mathbf{C}$). 2. Singular, non-positive
definite matrices where classical inverses fail.

The standard solution involves regularization techniques (like the Ridge
penalty `lambda` built into
[`gw_lmsi()`](https://zankrut20.github.io/selection.index/reference/gw_lmsi.md)),
ensuring matrix invertibility at the cost of slight coefficient bias.

## Conclusion

The `selection.index` package extends base index selection directly into
the molecular realm. The **LMSI**
([`lmsi()`](https://zankrut20.github.io/selection.index/reference/lmsi.md))
allows breeders to leverage an established pre-calculated marker
scorings regression model. The **GW-LMSI**
([`gw_lmsi()`](https://zankrut20.github.io/selection.index/reference/gw_lmsi.md))
embraces the Genomic Selection paradigm by swallowing all marker and
phenotypic data jointly in a single predictive framework, utilizing
ridge regularization to bypass classical $p \gg n$ multicollinearity
obstacles.

When implemented on traits exhibiting low heritability, these
marker-assisted pipelines can drastically outcompete standard phenotypic
indices in maximizing Expected Genetic Gain per selection cycle.
