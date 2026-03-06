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
  Justification: The installed size is primarily driven by the `libs` directory (~8.0Mb). This is expected and due to the heavy C++ template metaprogramming from `RcppEigen` used in the core mathematical engine.

## Reverse dependencies
I have checked 0 reverse dependencies and found no issues.

---

## Submission Summary
This is a resubmission to address the CRAN check failures reported in version 2.0.0 (requested to be fixed by 2026-03-16). 

**Fixes for CRAN:**
* **Fixed ERROR on macOS ARM64 (r-oldrel-macos-arm64):** Fixed the `invalid comparison with complex values` testing error in `ppg_esim`. The `eigen()` function returned tiny imaginary parts on macOS ARM64 for non-symmetric matrices; we now wrap the values in `Re()` before rank evaluation.

**Other Enhancements in 2.0.1:**
* **Expanded CI Coverage:** Integrated a comprehensive 17-runner GitHub Actions matrix mirroring all official CRAN check flavors including explicit ARM64 support for macOS, Windows, and Linux.
* **Documentation:** Automated spell-checking integrated into the CI pipeline via the `spelling` package.

Thank you to the CRAN team for reviewing this update.