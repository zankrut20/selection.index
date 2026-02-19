# Desired Gains Index (DG-LPSI)

Implements the Pesek & Baker (1969) Desired Gains Index where breeders
specify target genetic gains instead of economic weights. This enhanced
version includes calculation of implied economic weights and feasibility
checking.

## Usage

``` r
dg_lpsi(
  pmat,
  gmat,
  d,
  return_implied_weights = TRUE,
  check_feasibility = TRUE,
  selection_intensity = 2.063
)
```

## Arguments

- pmat:

  Phenotypic variance-covariance matrix (n_traits x n_traits)

- gmat:

  Genotypic variance-covariance matrix (n_traits x n_traits)

- d:

  Vector of desired genetic gains (length n_traits). Example: d = c(1.5,
  0.8, -0.2) means gain +1.5 in trait 1, +0.8 in trait 2, -0.2 in trait
  3.

- return_implied_weights:

  Logical - calculate implied economic weights? (default: TRUE)

- check_feasibility:

  Logical - warn if desired gains are unrealistic? (default: TRUE)

- selection_intensity:

  Selection intensity i (default: 2.063)

## Value

List with:

- `summary` - Data frame with coefficients, metrics, and implied weights

- `b` - Vector of selection index coefficients

- `Delta_G` - Named vector of achieved genetic gains per trait

- `desired_gains` - Named vector of desired gains (input d)

- `gain_errors` - Difference between desired and achieved gains

- `implied_weights` - Economic weights that would achieve these gains in
  Smith-Hazel LPSI

- `implied_weights_normalized` - Normalized implied weights (max
  absolute = 1)

- `feasibility` - Data frame with feasibility analysis per trait

- `hI2` - Index heritability

- `rHI` - Index accuracy

## Details

**Mathematical Formulation:**

1\. Index coefficients: \\\mathbf{b} = \mathbf{G}^{-1}\mathbf{d}\\

2\. Expected response: \\\Delta \mathbf{G} = (i/\sigma_I)
\mathbf{G}\mathbf{b}\\

**CRITICAL: Scale Invariance Property**

The achieved gains \\\Delta\mathbf{G}\\ are determined by selection
intensity (i), genetic variance (G), and phenotypic variance (P), NOT by
scaling \\\mathbf{b}\\. If you multiply \\\mathbf{b}\\ by constant c,
\\\sigma_I\\ also scales by c, causing complete cancellation in
\\\Delta\mathbf{G} = (i/(c\sigma_I))\mathbf{G}(c\mathbf{b}) =
(i/\sigma_I)\mathbf{G}\mathbf{b}\\.

**What DG-LPSI Actually Achieves:**

\- Proportional gains matching the RATIOS in d (not absolute
magnitudes) - Achieved magnitude depends on biological/genetic
constraints - Use feasibility checking to verify if desired gains are
realistic

3\. Implied economic weights (Section 1.4 of Chapter 4):
\$\$\hat{\mathbf{w}} = \mathbf{G}^{-1} \mathbf{P} \mathbf{b}\$\$

The implied weights represent the economic values that would have been
needed in a Smith-Hazel index to achieve the desired gain PROPORTIONS.
Large implied weights indicate traits that are "expensive" to improve
(low heritability or unfavorable correlations), while small weights
indicate traits that are "cheap" to improve.

**Feasibility Checking:**

The function estimates maximum possible gains as approximately 3.0 \*
sqrt(G_ii) (assuming very intense selection with i ~ 3.0) and warns if
desired gains exceed 80% of these theoretical maxima.

## References

Pesek, J., & Baker, R. J. (1969). Desired improvement in relation to
selection indices. *Canadian Journal of Plant Science*, 49(6), 803-804.

## Examples

``` r
# Load data
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Specify desired gains (e.g., increase each trait by 1 unit)
desired_gains <- rep(1, ncol(gmat))

# Calculate Desired Gains Index with all enhancements
result <- dg_lpsi(pmat, gmat, d = desired_gains)
#> Warning: Desired gains may be unrealistic for trait(s): rpp, ppr, spp, pw
#> Desired gains exceed 80% of theoretical maximum (i = 2.06).
#> Consider reducing targets or using higher selection intensity.

# View summary
print(result$summary)
#>       b.1    b.2     b.3      b.4    b.5       b.6     b.7 Delta_G     hI2 rHI
#> 1 11.4784 3.8899 -17.521 -17.3052 0.5096 -115.2781 71.4102 53.5871 -0.0931  NA
#>   implied_w.1 implied_w.2 implied_w.3 implied_w.4 implied_w.5 implied_w.6
#> 1    -46.6372     14.1965    178.2319     72.9888    -26.4785   -256.3287
#>   implied_w.7 implied_w_norm.1 implied_w_norm.2 implied_w_norm.3
#> 1    738.7446          -0.0631           0.0192           0.2413
#>   implied_w_norm.4 implied_w_norm.5 implied_w_norm.6 implied_w_norm.7
#> 1           0.0988          -0.0358           -0.347                1

# Extract implied weights to understand relative "cost" of gains
print(result$implied_weights_normalized)
#>        sypp         dtf         rpp         ppr         ppp         spp 
#> -0.06313034  0.01921701  0.24126326  0.09880111 -0.03584256 -0.34697876 
#>          pw 
#>  1.00000000 

# Check feasibility
print(result$feasibility)
#>      trait desired_gain achieved_gain genetic_sd max_possible_gain
#> sypp  sypp            1    0.07942154  1.1209826            2.3126
#> dtf    dtf            1    0.07942154  1.2490710            2.5768
#> rpp    rpp            1    0.07942154  0.3639840            0.7509
#> ppr    ppr            1    0.07942154  0.4931853            1.0174
#> ppp    ppp            1    0.07942154  0.9801869            2.0221
#> spp    spp            1    0.07942154  0.1319506            0.2722
#> pw      pw            1    0.07942154  0.1015121            0.2094
#>      feasibility_ratio is_realistic
#> sypp            0.4324         TRUE
#> dtf             0.3881         TRUE
#> rpp             1.3317        FALSE
#> ppr             0.9829        FALSE
#> ppp             0.4945         TRUE
#> spp             3.6736        FALSE
#> pw              4.7751        FALSE
```
