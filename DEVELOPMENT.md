# Quick Development Guide for selection.index

## First-Time Setup

``` r
# 1. Install all development dependencies
source("dev-setup.R")

# 2. Load the package
devtools::load_all()

# 3. Run tests to verify everything works
devtools::test()
```

## Daily Development Workflow

``` r
# Load your changes
devtools::load_all()

# Write code...

# Update documentation
devtools::document()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Common Commands

### Package Loading & Building

``` r
devtools::load_all()           # Load package for interactive use
devtools::install()            # Install package locally
devtools::build()              # Build package tarball
```

### Documentation

``` r
devtools::document()           # Generate .Rd files from roxygen comments
pkgdown::build_site()          # Build package website
rmarkdown::render("README.Rmd") # Update README.md
```

### Testing

``` r
devtools::test()               # Run all tests
devtools::test_coverage()      # Check test coverage
testthat::test_file("tests/testthat/test-*.R") # Run specific test
```

### Quality Checks

``` r
devtools::check()              # Comprehensive R CMD check
lintr::lint_package()          # Code linting
styler::style_pkg()            # Format code
goodpractice::gp()             # Best practices check
```

### CRAN Preparation

``` r
devtools::check()              # Local check
rhub::check_for_cran()         # Multi-platform check
devtools::check_win_devel()    # Windows devel check
devtools::check_win_release()  # Windows release check
```

## File Locations

- **R Code:** `R/*.R`
- **Tests:** `tests/testthat/test-*.R`
- **Documentation:** `man/*.Rd` (auto-generated)
- **Vignettes:** `vignettes/*.Rmd`
- **Data:** `data/*.rda`, `data-raw/*.R`

## Need Help?

- See `DEPENDENCIES.md` for detailed dependency information
- Run `?usethis` for package development helpers
- Visit <https://r-pkgs.org/> for comprehensive guidance
