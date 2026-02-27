# Changelog

## selection.index 2.0.0

### Breaking changes

- All primary user-facing functions have transitioned from dot-notation
  to snake_case notation for better clarity and to avoid S3 method
  dispatch conflicts.
  - `phen.varcov()` is now
    [`phen_varcov()`](https://zankrut20.github.io/selection.index/reference/phen_varcov.md)
  - `gen.varcov()` is now
    [`gen_varcov()`](https://zankrut20.github.io/selection.index/reference/gen_varcov.md)
  - `weight.mat()` is now
    [`weight_mat()`](https://zankrut20.github.io/selection.index/reference/weight_mat.md)
  - `gen.advance()` is now
    [`gen_advance()`](https://zankrut20.github.io/selection.index/reference/gen_advance.md)
- The base combinatorial selection functions have been replaced with the
  dedicated
  [`lpsi()`](https://zankrut20.github.io/selection.index/reference/lpsi.md)
  function (Linear Phenotypic Selection Index).

### Major changes

- The variance-covariance calculation engine, system solvers, and
  combinatorial index builders are now fully powered by `Rcpp` and
  `RcppEigen`.
- Index evaluations for large trait configurations are now exponentially
  faster due to primitive matrix operations.

### New features

- Introduced new suites for evaluating linear genomic selection indices
  ([`lgsi()`](https://zankrut20.github.io/selection.index/reference/lgsi.md))
  and combining genomic/phenotypic data
  ([`ppg_lgsi()`](https://zankrut20.github.io/selection.index/reference/ppg_lgsi.md)).
- Introduced new modules for tracking continuous index performance
  across multiple stages of breeding trials
  ([`mlpsi()`](https://zankrut20.github.io/selection.index/reference/mlpsi.md),
  [`mlgsi()`](https://zankrut20.github.io/selection.index/reference/mlgsi.md)).
- Implemented specialized solvers to maximize genetic advance while
  specifically constraining genetic gain for restrictive/undesired
  traits to zero
  ([`rlpsi()`](https://zankrut20.github.io/selection.index/reference/rlpsi.md),
  [`dg_lpsi()`](https://zankrut20.github.io/selection.index/reference/dg_lpsi.md)).
- Adopted Eigen-decomposition based selection methods for phenotypic
  ([`esim()`](https://zankrut20.github.io/selection.index/reference/esim.md))
  and genomic
  ([`gesim()`](https://zankrut20.github.io/selection.index/reference/gesim.md))
  evaluations.
- Added specialized infrastructure for evaluating linear marker
  selection indices
  ([`lmsi()`](https://zankrut20.github.io/selection.index/reference/lmsi.md)).
- Added new comprehensive toolsets to actively simulate and visualize
  multi-cycle genetic advance over time under varied selection index
  intensities and environmental variances
  ([`simulate_selection_cycles()`](https://zankrut20.github.io/selection.index/reference/simulate_selection_cycles.md)).
- Created a complete `pkgdown` official website.
- Added a complete vignette suite detailing the mathematical foundation
  and code usage for all new marker, genomic, multi-stage, constrained,
  and stochastic features.

## selection.index 1.2.1

CRAN release: 2026-01-26

- improve the overall performance of the package

## selection.index 1.2.0

CRAN release: 2023-09-19

- add new function `meanPerformance()` for calculating the mean
  performance in randomized block design data
- remove some typos from the documentation

## selection.index 1.1.4

CRAN release: 2022-06-13

- removed two functions `sel.index()` and `sel.score.rank()` for easy
  implementation of `comb.indices()`

## selection.index 1.1.3

CRAN release: 2021-09-25

- Added a new function `gen.advance()` for genetic advance calculation

## selection.index 1.1.2

CRAN release: 2021-07-26

- rank column added in `rcomb.indices()` for ranking of index

## selection.index 1.1.1

CRAN release: 2021-06-03

- Added a new function `rcomb.indices()` for calculated possible
  selection indices excluding single character or particular index from
  all possible selection indices

## selection.index 1.1.0

CRAN release: 2021-04-23

- Bug fixes in `phen.varcov()`
- Added a new function `comb.indices()` for calculating possible
  selection indices with the group/pairs of traits/characters
- removed function: `rank.index()`

## selection.index 1.0.0

CRAN release: 2021-04-06

- Added a `NEWS.md` file to track changes to the package.
