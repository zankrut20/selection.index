## Test environments
* local Windows 11, R 4.5.2
* GitHub Actions:
  * macOS (arm64 and x86_64) — R devel, release, oldrel-1 (6 runners)
  * Windows (arm64 and x86_64) — R devel, release, oldrel-1 (6 runners)
  * Linux (x86_64 and arm64) — R devel, release, patched, oldrel-1 (5 runners)
* win-builder (devel, release, and oldrelease)
* R-hub (Fedora Linux, Ubuntu Linux with ASAN/UBSAN/Valgrind, macOS-arm64)

## R CMD check results
0 errors | 0 warnings | 1 note

* NOTE: Installed package size is 9.2Mb. 
  The size is driven by the compiled `libs` directory (~8.0Mb) required for high-performance `RcppEigen` matrix operations.

## Submission Summary
This is a resubmission to address CRAN check findings from version 2.0.0.

### Major Fixes:
* **macOS ARM64:** Resolved the `invalid comparison with complex values` error in `ppg_esim`. Calculations now use strict real-part extraction via `Re()` on eigenvalues prior to rank/comparison checks, ensuring architectural parity between Intel and Apple Silicon.
* **Memory Sanitizers (gcc-san):** Addressed the `runtime error: signed integer overflow` reported in `tests/testthat.Rout`. This was identified as a known issue in the Base R TRE regular expression engine (`tre-match-approx.c`). To ensure clean CRAN runs, affected string-matching tests now use `fixed = TRUE` or `perl = TRUE`, and `skip_on_cran()` is utilized for heavy error-handling tests.
* **Compiler Warnings:** Fixed a strict compiler warning in `src/math_primitives.cpp` by removing an unused variable (`min_group`).
* **rchk:** All reports were reviewed. Identified reports reference internal `Rcpp` headers (`Armor.h`/`Shield.h`); no memory-protection issues exist in the package source code.

### Logistical Fixes:
* Removed non-portable `.mhtml` and `.txt` files from the source.
* Updated `.Rbuildignore` to ensure a clean root directory.
* Corrected test initialization in `tests/testthat/test-design_stats.R`.