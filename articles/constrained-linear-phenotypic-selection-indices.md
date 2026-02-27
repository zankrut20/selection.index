# Constrained Linear Phenotypic Selection Indices

## Introduction

The linear phenotypic selection index (LPSI), the null restricted LPSI
(RLPSI), and the predetermined proportional gains LPSI (PPG-LPSI) are
the main phenotypic selection indices used to predict the net genetic
merit and select parents for the next selection cycle. The LPSI is an
unrestricted index, whereas the RLPSI and the PPG-LPSI allow
restrictions equal to zero and predetermined proportional gain
restrictions respectively to be imposed on the expected genetic gain
values of the trait.

One additional restricted index is the desired gains LPSI (DG-LPSI),
which does not require economic weights and, in a similar manner to the
PPG-LPSI, allows restrictions to be imposed on the expected genetic gain
values of the trait based on a predetermined level.

In this vignette, we demonstrate the theory and practical application of
RLPSI, PPG-LPSI, and DG-LPSI using the `selection.index` package.

### Sample Data Preparation

We use the synthetic `maize_pheno` dataset built into the package to
calculate the phenotypic and genotypic variance-covariance matrices.

``` r
library(selection.index)

# Load the built-in maize phenotype dataset
data("maize_pheno")

# Extract the traits of interest
traits <- c("Yield", "PlantHeight", "DaysToMaturity")

# Calculate Genotypic (gmat) and Phenotypic (pmat) covariance matrices
gmat <- gen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)
pmat <- phen_varcov(maize_pheno[, traits], maize_pheno$Genotype, maize_pheno$Block)

# Define economic weights for the three traits
wmat <- weight_mat(data.frame(Trait = traits, Weight = c(1, -1, 1)))
```

## 1. The Null Restricted Linear Phenotypic Selection Index (RLPSI)

The main objective of the RLPSI is to optimize, under some null
restrictions, the selection response, to predict the net genetic merit
$H = \mathbf{w}\prime\mathbf{g}$, and to select the individuals with the
highest net genetic merit values as parents of the next generation. The
RLPSI allows restrictions equal to zero to be imposed on the expected
genetic gains of some traits.

Assuming we want to minimize the mean squared difference between the
index $I$ and merit $H$ ($E\left\lbrack (H - I)^{2} \right\rbrack$)
under the restriction $\mathbf{C}\prime\mathbf{b} = \mathbf{0}$, we
minimize the following function:

$$\Psi(\mathbf{b},\mathbf{v}) = \mathbf{b}^{\prime}{\mathbf{P}\mathbf{b}} + \mathbf{w}^{\prime}{\mathbf{G}\mathbf{w}} - 2\mathbf{w}^{\prime}{\mathbf{G}\mathbf{b}} + 2\mathbf{v}^{\prime}\mathbf{C}^{\prime}\mathbf{b}$$

where $\mathbf{C}\prime = \mathbf{U}\prime\mathbf{G}$,
$\mathbf{U}\prime$ is the matrix of null restrictions, and $\mathbf{v}$
is a vector of Lagrange multipliers. The solution yields the RLPSI
vector of coefficients:

$$\mathbf{b}_{R} = {\mathbf{K}\mathbf{b}}$$

where $\mathbf{K} = \lbrack\mathbf{I} - \mathbf{Q}\rbrack$,
$\mathbf{Q} = \mathbf{P}^{- 1}\mathbf{C}\left( \mathbf{C}\prime\mathbf{P}^{- 1}\mathbf{C} \right)^{- 1}\mathbf{C}\prime$,
and $\mathbf{b} = \mathbf{P}^{- 1}{\mathbf{G}\mathbf{w}}$ is the LPSI
vector of coefficients.

### Practical Application: RLPSI

Using the `rlpsi` function, we can constrain the genetic gain of certain
traits to be strictly zero. For instance, we may want to maximize
overall gain while ensuring that Plant Height (trait \#2) does not
change.

``` r
# Restrict trait 2 (PHT) to have ZERO expected genetic gain
rlpsi_res <- rlpsi(
  pmat = pmat,
  gmat = gmat,
  wmat = wmat,
  restricted_traits = c(2)
)

# View the summary and coefficients
print(rlpsi_res$summary)
#>      b.1     b.2    b.3       GA      PRE  Delta_G    rHI    hI2
#> 1 0.0551 -3.3075 2.1743 108.2417 10824.17 108.2417 0.3153 0.0994

# View the expected genetic gains (Delta_G)
print(rlpsi_res$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>   1.059783e+02  -1.827973e-08   2.263426e+00
```

Notice that the expected genetic gain (`Delta_G`) for the 2nd trait is
effectively zero.

## 2. The Predetermined Proportional Gains Linear Phenotypic Selection Index (PPG-LPSI)

Unlike the RLPSI, the predetermined proportional gains phenotypic
selection index (PPG-LPSI) allows restrictions *different from zero* to
be imposed, ensuring that traits gain in strictly predefined proportions
to each other.

Minimizing $E\left\lbrack (H - I)^{2} \right\rbrack$ under the
proportional restriction $\mathbf{M}\prime\mathbf{b} = 0$, where
$\mathbf{M}\prime = \mathbf{D}\prime\mathbf{C}\prime$ and
$\mathbf{D}\prime$ is a Mallard (1972) matrix, leads to:

$$\Phi(\mathbf{b},\mathbf{v}) = \mathbf{b}^{\prime}{\mathbf{P}\mathbf{b}} + \mathbf{w}^{\prime}{\mathbf{G}\mathbf{w}} - 2\mathbf{w}^{\prime}{\mathbf{G}\mathbf{b}} + 2\mathbf{v}^{\prime}\mathbf{M}^{\prime}\mathbf{b}$$

The vector that minimizes the mean squared error under this restriction
is:

$$\mathbf{b}_{M} = \mathbf{K}_{M}\mathbf{b}$$

where
$\mathbf{K}_{M} = \left\lbrack \mathbf{I} - \mathbf{Q}_{M} \right\rbrack$
and
$\mathbf{Q}_{M} = \mathbf{P}^{- 1}\mathbf{M}\left( \mathbf{M}\prime\mathbf{P}^{- 1}\mathbf{M} \right)^{- 1}\mathbf{M}\prime$.
Alternatively, Tallis (1985) formulated it using a proportionality
constant $\theta$:

$$\mathbf{b}_{T} = \mathbf{b}_{R} + \theta\mathbf{δ}$$

### Practical Application: PPG-LPSI

The `ppg_lpsi` function enforces these desired proportional gains.
Suppose we want the gains across the three traits to follow the
proportion `2 : 1 : 1`.

``` r
# Specify the desired proportions
k_proportions <- c(2, 1, 1)

# Calculate the PPG-LPSI
ppg_res <- ppg_lpsi(pmat = pmat, gmat = gmat, k = k_proportions, wmat = wmat)

# View the expected genetic gains
print(ppg_res$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>      700636.74       18363.40       10206.29
```

If we observe the resulting `Delta_G` values, their relative proportions
should approximate the `2 : 1 : 1` ratio dictated by `k`.

## 3. The Desired Gains Linear Phenotypic Selection Index (DG-LPSI)

The desired gains linear phenotypic selection index (DG-LPSI) is unique
in that it *does not require economic weights* $\mathbf{w}$. Instead,
the breeder specifies the exact desired target genetic gains
$\mathbf{d}$ directly.

Because the expected genetic gain is
$\mathbf{E} = k_{I}\frac{\mathbf{G}\mathbf{b}}{\sigma_{I}}$, if we set
${\mathbf{G}\mathbf{b}} = \mathbf{d}$, then we want to minimize the
variance of the index $\sigma_{I}$ subject to
${\mathbf{G}\mathbf{b}} = \mathbf{d}$ to maximize $\mathbf{E}$:

$$\Phi_{DG}(\mathbf{b},\mathbf{v}) = 0.5\left( \mathbf{b}^{\prime}{\mathbf{P}\mathbf{b}} \right) + \mathbf{v}^{\prime}\left( {\mathbf{G}\mathbf{b}} - \mathbf{d} \right)$$

Solving this yields the DG-LPSI vector of coefficients:

$$\mathbf{b}_{DG} = \mathbf{P}^{- 1}\mathbf{G}\left( {\mathbf{G}\mathbf{P}}^{- 1}\mathbf{G} \right)^{- 1}\mathbf{d}$$

### Practical Application: DG-LPSI

We can run the `dg_lpsi` function solely based on desired absolute
gains. Let’s aim to increase Yield by 5 units, decrease Plant Height by
2 units, and increase Days to Maturity by 1 unit.

``` r
# Explicit vector of desired absolute genetic gains
desired_gains <- c(5, -2, 1)

# Calculate DG-LPSI
dg_res <- dg_lpsi(pmat = pmat, gmat = gmat, d = desired_gains)

# Check the achieved proportional genetic gains
print(dg_res$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>      3.3431312     -1.3372525      0.6686262

# The DG-LPSI also calculates implied Smith-Hazel economic weights
print(dg_res$implied_weights_normalized)
#>          Yield    PlantHeight DaysToMaturity 
#>    -0.01969999     0.19583750     1.00000000
```

The output gives us the selection index coefficients $\mathbf{b}$ that
would achieve the target gains, as well as the implied economic weights
showing the “absolute cost” or relative importance required to achieve
those gains.

## 4. Advanced Concepts and Trade-offs

### 4.1 Statistical Properties and Efficiency Trade-offs

Imposing restrictions on an index inherently limits its potential to
maximize the overall selection response $R_{R}$. The `selection.index`
package allows you to evaluate this “cost of restriction” by calculating
both the relative selection efficiency (`PRE` - Percentage of Response
Efficiency) and the heritability (`hI2`) of the constrained index versus
the unconstrained equivalent.

When reviewing the `summary` data frame output of our indices above: -
**`rHI`**: The correlation between the constrained index and the true
net genetic merit ($H$). As more restrictions are added, this
correlation decreases. - **`PRE`**: Shows the expected genetic advance
relative to a base index. Constrained indices will universally show
lower `PRE` than their unconstrained LPSI counterpart.

### 4.2 Manual Constraint Matrix Construction

By default,
[`rlpsi()`](https://zankrut20.github.io/selection.index/reference/rlpsi.md)
and
[`ppg_lpsi()`](https://zankrut20.github.io/selection.index/reference/ppg_lpsi.md)
accept user-friendly arguments like `restricted_traits = c(2)` or
`k = c(2, 1, 1)` to automatically construct the mathematical constraint
matrices behind the scenes.

However, advanced users may wish to manually define these matrices,
particularly for complex scenarios.

For the **RLPSI**, the constraint matrix $\mathbf{C}$ is a $t \times r$
matrix (where $t$ is the number of traits and $r$ is the number of
restrictions). Each column represents a restricted trait with a `1` at
the restricted index and `0`s elsewhere:

``` r
# Manually restrict traits #1 and #3:
# Create a 3x2 constraint matrix
C_matrix <- matrix(
  c(
    1, 0, 0, # Restrict trait 1
    0, 0, 1
  ), # Restrict trait 3
  nrow = 3, ncol = 2
)

# Pass directly to RLPSI
rlpsi_manual <- rlpsi(pmat = pmat, gmat = gmat, wmat = wmat, C = C_matrix)
print(rlpsi_manual$Delta_G)
#>          Yield    PlantHeight DaysToMaturity 
#>    -0.01072552    -2.76848255     0.74035592
```

The
[`ppg_lpsi()`](https://zankrut20.github.io/selection.index/reference/ppg_lpsi.md)
function automatically manages the complex proportional constraint
matrix based on your inputs $\mathbf{k}$, guaranteeing the proportional
relationship $\Delta\mathbf{G} = \phi\mathbf{k}$.

### 4.3 Applicability Limitations (Hazel, 1943)

While these indices are powerful predictive models, their real-world
applicability operates under the limitations established by classical
breeding theory (Hazel, 1943):

1.  **Economic Variability:** Relative economic values ($\mathbf{w}$)
    fluctuate by location, market demand, and enterprise goals.
2.  **Genetic Context:** A given genetic covariance matrix
    ($\mathbf{G}$) is specific to the mating system and breeding
    population evaluated.
3.  **Environmental Variance:** Phenotypic standard deviations
    ($\mathbf{P}$) differ among environments and management practices.
4.  **Sampling Error:** Small datasets fail to provide the reliable
    genetic variance constants required to form precise predictive
    indices.

## Conclusion

The `selection.index` package provides highly flexible options beyond
the simple unrestricted index. The **RLPSI** is useful when strictly
zero-change bounds are needed on specific traits. The **PPG-LPSI**
caters to breeders focusing on a strict multi-trait improvement ratio.
Finally, the **DG-LPSI** provides a powerful alternative for scenarios
where economic weights are notoriously difficult to estimate, but ideal
target gains are known. When using these constrained models, breeders
must continuously balance their specific phenotypic goals against the
intrinsic statistical costs to overall genetic efficiency.
