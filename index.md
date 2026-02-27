# selection.index

Stop wasting days manually assembling complex variance-covariance
matrices and resolving selection weights. **`selection.index`** is the
definitive, production-ready R package built specifically for
agricultural statisticians and quantitative geneticists to instantly
calculate optimal multi-trait selection indices. Whether you are running
classical phenotypic field trials or building advanced genomic
evaluations, this package provides a fully integrated mathematical
engine that directly maximizes genetic advance and substantially
accelerates your breeding pipeline.

**📚 [Read the Official Documentation & Vignettes on our pkgdown
Website](https://zankrut20.github.io/selection.index/)**

## Installation

Install the stable release from CRAN:

``` r
install.packages("selection.index")
```

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("zankrut20/selection.index")
```

## Quick Start: Optimizing Phenotypic Selection

The package includes several pre-configured structural datasets
(`seldata`, `maize_pheno`, `maize_geno`) to test your models. Here, we
demonstrate a fundamental phenotypic index evaluating multiple field
traits simultaneously.

Before combining traits into a unified index, breeders must cleanly
partition the total phenotypic variation. We first use the
variance-covariance functions
([`phen_varcov()`](https://zankrut20.github.io/selection.index/reference/phen_varcov.md)
and
[`gen_varcov()`](https://zankrut20.github.io/selection.index/reference/gen_varcov.md))
to extract the precise phenotypic and genotypic matrices directly from
your replicated trial data.

``` r
# Using the built-in 'seldata' dataset (Replication in Col 1, Genotype in Col 2)
p_mat <- phen_varcov(data = seldata[, 3:9], genotypes = seldata[, 2], replication = seldata[, 1])
g_mat <- gen_varcov(data = seldata[, 3:9], genotypes = seldata[, 2], replication = seldata[, 1])
```

Next, we establish strict economic objectives for our evaluated traits.
The
[`weight_mat()`](https://zankrut20.github.io/selection.index/reference/weight_mat.md)
function rapidly converts raw economic weight arrays into the formalized
mathematical format required for the core index calculation.

``` r
# Convert the built-in 'weight' configuration matrix
w_mat <- weight_mat(weight)
```

With our covariance structures and economic priorities firmly defined,
we must evaluate the baseline potential of our population. We calculate
the expected genetic gain for our primary target trait (e.g., yield)
using the
[`gen_advance()`](https://zankrut20.github.io/selection.index/reference/gen_advance.md)
function.

``` r
# Calculate expected genetic advance for the primary trait
base_advance <- gen_advance(phen_mat = p_mat[1, 1], gen_mat = g_mat[1, 1], weight_mat = w_mat[1, 1])
base_advance
#>        [,1]
#> [1,] 1.2009
```

Finally, to make definitive selection decisions, we dynamically evaluate
trait combinations. The
[`lpsi()`](https://zankrut20.github.io/selection.index/reference/lpsi.md)
function constructs robust combinatorial selection indices and
immediately ranks them by predicted relative efficiency, allowing you to
confidently select the mathematically superior combination.

``` r
# Calculate and rank combinatorial selection indices
optimal_index <- lpsi(
  ncomb = 1, pmat = p_mat, gmat = g_mat,
  wmat = weight[, 2:3], wcol = 1, GAY = base_advance
)
head(optimal_index)
#>   ID    b.1     GA      PRE Delta_G    rHI    hI2 Rank
#> 1  1 0.6316 2.8125 234.2021  2.8125 0.4825 0.2697    1
#> 2  2 0.2396 1.3036 108.5481  1.3036 0.2237 0.2242    4
#> 3  3 1.4869 2.1526 179.2462  2.1526 0.3693 0.2690    2
#> 4  4 0.4507 0.9084  75.6402  0.9084 0.1558 0.2549    6
#> 5  5 0.3164 1.6236 135.1988  1.6236 0.2786 0.1553    3
#> 6  6 1.4499 1.0292  85.7009  1.0292 0.1766 0.1471    5
```

*(For advanced modeling involving marker-assisted indices, multistage
selection, or genomic predictions, please load the included
`maize_pheno` and `maize_geno` datasets and consult our [pkgdown
vignettes](https://zankrut20.github.io/selection.index/).)*

## Citation

We rely directly on academic citations to maintain and justify the
continued development of this software ecosystem. If you use
`selection.index` in a publication, breeding pipeline, or research
project, it is imperative that you cite it:

``` r
citation("selection.index")
```

## Support & Consulting

### Direct-to-Client Consulting

Are you scaling up a commercial breeding operation or standardizing a
massive academic quantitative pipeline? I offer dedicated, expert
consulting to researchers tailored for custom selection index
development, stochastic simulation modeling, and full software pipeline
integration.

**Note: I operate strictly direct-to-client. There are absolutely no
middlemen, generic agencies, or intermediaries involved. You work
directly with the package creator.**

To discuss consulting scope and availability, contact:
`zankrut20@gmail.com`

### Tip Jar

If this package has saved you substantial time, simplified your data
pipeline, or helped you publish faster, please consider supporting its
independent development! Monetary contributions directly fund repository
maintenance and the addition of modern quantitative genetic models.

- **Indian Users (UPI):** `zankrut20@upi` *(Replace with proper UPI)*
- **International Users:** Available Soon
