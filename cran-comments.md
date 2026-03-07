## Test environments
* local Windows 11, R 4.5.2
* GitHub Actions:
  * macOS (arm64 and x86_64) — R devel, release, oldrel-1 (total 6 runners)
  * Windows (arm64 and x86_64) — R devel, release, oldrel-1 (total 6 runners)
  * Linux (x86_64 and arm64) — R devel, release, patched, oldrel-1 (total 5 runners)
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 1 note

* NOTE: Installed package size is 9.2Mb.
  The installed size is primarily driven by the `libs` directory (~8.0Mb) due to `RcppEigen`-based C++ code.

## Submission Summary
This is a resubmission to address CRAN check findings from version 2.0.0.

Changes made:
* Fixed macOS ARM64 test error (`invalid comparison with complex values`) in `ppg_esim` by using `Re()` on eigenvalues before rank checks.
* Removed non-portable `.mhtml` files and updated `.Rbuildignore` to exclude development log files (`.txt`, `.mhtml`).
* Updated string-matching tests to use explicit matching (`fixed = TRUE` or `perl = TRUE`) and kept `skip_on_cran()` in CRAN-sensitive heavy/error-handling tests.
* Removed unused C++ variable in `src/math_primitives.cpp` (`min_group`) to resolve strict compiler warning checks.
* Corrected test setup in `tests/testthat/test-design_stats.R` for invalid `design_type` test input initialization.

rchk note:
* Reports referencing `Rcpp` internal headers (`Armor.h` / `Shield.h`) were reviewed; no corresponding memory-protection issue was identified in package source files.