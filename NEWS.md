# selection.index (development version)

# selection.index 2.0.0

## Breaking changes

* All primary user-facing functions have transitioned from dot-notation to snake_case notation for better clarity and to avoid S3 method dispatch conflicts.
  * `phen.varcov()` is now `phen_varcov()`
  * `gen.varcov()` is now `gen_varcov()`
  * `weight.mat()` is now `weight_mat()`
  * `gen.advance()` is now `gen_advance()`
* The base combinatorial selection functions have been replaced with the dedicated `lpsi()` function (Linear Phenotypic Selection Index).

## Major changes

* The variance-covariance calculation engine, system solvers, and combinatorial index builders are now fully powered by `Rcpp` and `RcppEigen`.
* Index evaluations for large trait configurations are now exponentially faster due to primitive matrix operations.

## New features

* Introduced new suites for evaluating linear genomic selection indices (`lgsi()`) and combining genomic/phenotypic data (`ppg_lgsi()`).
* Introduced new modules for tracking continuous index performance across multiple stages of breeding trials (`mlpsi()`, `mlgsi()`).
* Implemented specialized solvers to maximize genetic advance while specifically constraining genetic gain for restrictive/undesired traits to zero (`rlpsi()`, `dg_lpsi()`).
* Adopted Eigen-decomposition based selection methods for phenotypic (`esim()`) and genomic (`gesim()`) evaluations.
* Added specialized infrastructure for evaluating linear marker selection indices (`lmsi()`).
* Added new comprehensive toolsets to actively simulate and visualize multi-cycle genetic advance over time under varied selection index intensities and environmental variances (`simulate_selection_cycles()`).
* Created a complete `pkgdown` official website.
* Added a complete vignette suite detailing the mathematical foundation and code usage for all new marker, genomic, multi-stage, constrained, and stochastic features.

# selection.index 1.2.1

* improve the overall performance of the package

# selection.index 1.2.0

* add new function `meanPerformance()` for calculating the mean performance in randomized block design data 
* remove some typos from the documentation

# selection.index 1.1.4

* removed two functions `sel.index()` and `sel.score.rank()` for easy implementation of `comb.indices()`

# selection.index 1.1.3

* Added a new function `gen.advance()` for genetic advance calculation

# selection.index 1.1.2

* rank column added in `rcomb.indices()` for ranking of index

# selection.index 1.1.1

* Added a new function `rcomb.indices()` for calculated possible selection indices excluding single character or particular index from all possible selection indices

# selection.index 1.1.0

* Bug fixes in `phen.varcov()` 
* Added a new function `comb.indices()` for calculating possible selection indices with the group/pairs of traits/characters
* removed function: `rank.index()` 

# selection.index 1.0.0

* Added a `NEWS.md` file to track changes to the package.
