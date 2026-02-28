## Test environments
* OS X (on GitHub Actions), R 4.3.2
* ubuntu 20.04 (on GitHub Actions), R 4.3.2
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Reverse dependencies
I have checked 0 reverse dependencies and found no issues.

---

## Submission Summary
This is a major release (1.2.1 -> 2.0.0) reflecting a comprehensive rewrite of the mathematical engine, sweeping new feature additions for modern quantitative genetics, and API rationalization.

**Major Changes in 2.0.0:**
* **C++ Integration:** The core mathematical solvers (matrices, combinatorial index permutations) have been transitioned to C++ using `Rcpp` and `RcppEigen` for a structural performance overhaul.
* **API Standardization:** Core functions transitioned from dot.notation to snake_case (e.g. `phen.varcov()` -> `phen_varcov()`) to comply with CRAN S3 dispatch recommendations and modernize the API.
* **Modern Genetics Methodologies Added:** Substantial new functions handling Genomic Selection Indices, Marker-Assisted Selection, Eigen-based indices, Constrained/Restricted indices, and Multi-Stage selection.
* **Simulation:** Added stochastic genetic advance simulation models explicitly tracking progressive selection over multiple generations.
* **Vignettes:** 10 comprehensive vignettes added mathematically documenting the index formulations.

Thank you to the CRAN team for reviewing this significant update.
