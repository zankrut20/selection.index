## Test environments
* local OS X install, R 4.3.2
* ubuntu 20.04 (on GitHub Actions), R 4.3.2
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Reverse dependencies
I have checked 0 reverse dependencies and found no issues.

---

## Submission Summary
This is a patch release (1.2.0 -> 1.2.1).

**Changes:**
* Major performance improvements to internal functions using optimized base R code.
* No new functionality or user-facing API changes were made.
* All external function signatures remain identical to the previous version.
