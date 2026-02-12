# Package Dependencies Documentation

This document provides a comprehensive overview of all dependencies used
in the `selection.index` package. This information is crucial for: -
Setting up a new development environment - Understanding what each
package is used for - Maintaining CRAN compliance - Future development
and maintenance

------------------------------------------------------------------------

## Table of Contents

1.  [Runtime Dependencies](#runtime-dependencies)
2.  [Development Dependencies](#development-dependencies)
3.  [Testing Dependencies](#testing-dependencies)
4.  [Documentation Dependencies](#documentation-dependencies)
5.  [Optional Quality Tools](#optional-quality-tools)
6.  [Quick Setup](#quick-setup)
7.  [CRAN Compliance Notes](#cran-compliance-notes)

------------------------------------------------------------------------

## Runtime Dependencies

These packages are required for the package to function and are listed
in the `Imports` field of `DESCRIPTION`. They will be automatically
installed when users install `selection.index`.

| Package | Version | Purpose                  | Functions Used                                                                                                                         |
|---------|---------|--------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| `stats` | Base R  | Statistical calculations | [`qt()`](https://rdrr.io/r/stats/TDist.html), [`pf()`](https://rdrr.io/r/stats/Fdist.html) for statistical tests in mean_performance.R |
| `utils` | Base R  | Utility functions        | [`combn()`](https://rdrr.io/r/utils/combn.html) for generating combinations in comb_indices.R                                          |

**Note:** Both `stats` and `utils` are base R packages that come
pre-installed with R, but they must be declared in DESCRIPTION for CRAN
compliance.

------------------------------------------------------------------------

## Development Dependencies

These packages are used during package development but are not required
by end users. They are listed in the `Suggests` field or installed
separately via `dev-setup.R`.

### Essential Development Tools

| Package     | Purpose                                   | When to Use                                        |
|-------------|-------------------------------------------|----------------------------------------------------|
| `devtools`  | Complete package development workflow     | Load package during development, run checks, build |
| `usethis`   | Automate package development tasks        | Create new files, functions, tests                 |
| `roxygen2`  | Generate documentation from code comments | Update .Rd files from @-tags in R code             |
| `pkgdown`   | Build package website                     | Create/update GitHub Pages site                    |
| `rcmdcheck` | Run R CMD check from R                    | Check package before submission                    |
| `remotes`   | Install packages from various sources     | Install dependencies from GitHub                   |

### Installation

``` r
install.packages(c("devtools", "usethis", "roxygen2", "pkgdown", "rcmdcheck", "remotes"))
```

------------------------------------------------------------------------

## Testing Dependencies

Required for running the package test suite.

| Package    | Version Required | Purpose                | Configuration              |
|------------|------------------|------------------------|----------------------------|
| `testthat` | \>= 3.0.0        | Unit testing framework | Tests in `tests/testthat/` |
| `covr`     | \-               | Code coverage analysis | Used in CI/CD pipelines    |

**Test Files:** - `test-comb_indices.R` - Tests for combination
indices - `test-design_stats.R` - Tests for design statistics -
`test-gen_advance.R` - Tests for genetic advance - `test-genphen.R` -
Tests for genotypic/phenotypic functions - `test-mean_performance.R` -
Tests for mean performance - `test-missing_value_estimation.R` - Tests
for missing value estimation - `test-rcomb_indices.R` - Tests for
reverse combination indices - `test-weight_mat.R` - Tests for weight
matrix

### Running Tests

``` r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-comb_indices.R")

# Check code coverage
covr::package_coverage()
```

------------------------------------------------------------------------

## Documentation Dependencies

Required for building vignettes and documentation.

| Package     | Purpose                        | Used In               |
|-------------|--------------------------------|-----------------------|
| `knitr`     | Dynamic report generation      | Vignettes, README.Rmd |
| `rmarkdown` | R Markdown document conversion | Vignettes, READMEs    |
| `markdown`  | Markdown rendering             | Documentation         |

**Documentation Files:** - `vignettes/Examples.Rmd` - Main package
vignette - `README.Rmd` - GitHub README source - `man/*.Rd` - Function
documentation (auto-generated by roxygen2)

### Building Documentation

``` r
# Update function documentation
devtools::document()

# Build vignettes
devtools::build_vignettes()

# Build pkgdown site
pkgdown::build_site()

# Compile README
rmarkdown::render("README.Rmd")
```

------------------------------------------------------------------------

## Optional Quality Tools

These packages help maintain code quality but are not required for basic
development.

| Package        | Purpose                   | Usage                   |
|----------------|---------------------------|-------------------------|
| `lintr`        | Static code analysis      | `lintr::lint_package()` |
| `styler`       | Automatic code formatting | `styler::style_pkg()`   |
| `goodpractice` | Best practices checker    | `goodpractice::gp()`    |

### Code Quality Checks

``` r
# Lint the package
lintr::lint_package()

# Format code consistently
styler::style_pkg()

# Check best practices
goodpractice::gp()
```

------------------------------------------------------------------------

## Quick Setup

### For New Developers

1.  **Clone the repository:**

    ``` bash
    git clone https://github.com/zankrut20/selection.index.git
    cd selection.index
    ```

2.  **Run the development setup script:**

    ``` r
    source("dev-setup.R")
    ```

3.  **Load the package for development:**

    ``` r
    devtools::load_all()
    ```

4.  **Run tests to ensure everything works:**

    ``` r
    devtools::test()
    ```

5.  **Check the package:**

    ``` r
    devtools::check()
    ```

### Manual Installation

If you prefer to install packages manually:

``` r
# Install all development dependencies
install.packages(c(
  # Development tools
  "devtools", "usethis", "roxygen2", "pkgdown", "rcmdcheck", "remotes",
  
  # Testing
  "testthat", "covr",
  
  # Documentation
  "knitr", "rmarkdown", "markdown",
  
  # Quality tools (optional)
  "lintr", "styler", "goodpractice"
))

# Install package dependencies
remotes::install_deps(dependencies = TRUE)
```

------------------------------------------------------------------------

## CRAN Compliance Notes

### Package Structure Compliance

1.  **DESCRIPTION File:**
    - ✅ All runtime dependencies in `Imports`
    - ✅ Development/testing dependencies in `Suggests`
    - ✅ R version requirement specified: `R (>= 2.10)`
    - ✅ Valid license: GPL (\>= 3)
    - ✅ All authors properly credited
    - ✅ Valid URLs and BugReports
2.  **Dependency Management:**
    - ✅ Base R packages (`stats`, `utils`) properly declared in Imports
    - ✅ No unnecessary dependencies
    - ✅ All `@importFrom` directives match DESCRIPTION
    - ✅ No entire packages imported (using selective imports)
3.  **Build Files:**
    - ✅ `.Rbuildignore` excludes development files
    - ✅ Development scripts not included in CRAN package
    - ✅ Documentation properly generated

### Dependency Best Practices

1.  **Minimal Runtime Dependencies:**
    - Only `stats` and `utils` (both base R packages)
    - No external CRAN dependencies for core functionality
    - Reduces maintenance burden and compatibility issues
2.  **Import Strategy:**
    - Use `@importFrom` for specific functions (not entire packages)
    - Clearly document what functions come from which packages
    - Keep NAMESPACE clean and minimal
3.  **Version Requirements:**
    - Specify minimum versions only when necessary
    - Current: `testthat (>= 3.0.0)` - for modern test features
    - Base R requirement: `R (>= 2.10)` - for LazyData support

### Pre-CRAN Submission Checklist

Before submitting to CRAN, run:

``` r
# 1. Update documentation
devtools::document()

# 2. Run comprehensive checks
devtools::check()

# 3. Check on R-hub (multiple platforms)
rhub::check_for_cran()

# 4. Check on Windows (if on Unix/Mac)
devtools::check_win_devel()
devtools::check_win_release()

# 5. Update NEWS.md with changes

# 6. Update cran-comments.md with test environments
```

------------------------------------------------------------------------

## Dependency Update Strategy

### When to Update Dependencies

1.  **Security vulnerabilities:** Update immediately
2.  **Bug fixes in dependencies:** Update if they affect functionality
3.  **New features needed:** Evaluate and update cautiously
4.  **Breaking changes:** Test thoroughly before updating

### How to Update

``` r
# Check for outdated packages
old.packages()

# Update all packages
update.packages(ask = FALSE)

# Or update selectively
install.packages("package_name")

# After updates, always run
devtools::test()
devtools::check()
```

### Monitoring Dependencies

- Review DESCRIPTION periodically
- Check for deprecated functions
- Monitor R package news and updates
- Test with latest R-devel periodically

------------------------------------------------------------------------

## Troubleshooting

### Common Issues

1.  **Missing dependencies during devtools::check():**

    ``` r
    remotes::install_deps(dependencies = TRUE)
    ```

2.  **Roxygen version mismatch:**

    ``` r
    install.packages("roxygen2")
    devtools::document()
    ```

3.  **Test failures after updates:**

    ``` r
    devtools::test()  # Check specific failures
    devtools::load_all()  # Reload package
    ```

4.  **Vignette build failures:**

    ``` r
    # Ensure all documentation packages installed
    install.packages(c("knitr", "rmarkdown", "markdown"))
    devtools::build_vignettes()
    ```

------------------------------------------------------------------------

## Resources

- [R Packages Book](https://r-pkgs.org/) - Comprehensive guide to
  package development
- [CRAN Repository
  Policy](https://cran.r-project.org/web/packages/policies.html)
- [Writing R
  Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
- [Package Development Cheat
  Sheet](https://github.com/rstudio/cheatsheets/blob/main/package-development.pdf)

------------------------------------------------------------------------

## Maintenance Log

| Date       | Action                   | Details                                        |
|------------|--------------------------|------------------------------------------------|
| 2026-02-05 | Initial Documentation    | Created comprehensive dependency documentation |
| 2026-02-05 | Added `stats` to Imports | Added missing base R package to DESCRIPTION    |
| 2026-02-05 | Created dev-setup.R      | Automated development environment setup        |

------------------------------------------------------------------------

**Last Updated:** February 5, 2026  
**Package Version:** 1.2.1.9000  
**Maintainer:** Zankrut Goyani <zankrut20@gmail.com>
