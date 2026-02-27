# =============================================================================
# Development Environment Setup for selection.index Package
# =============================================================================
#
# This script installs all packages required for developing, testing, and
# building the selection.index package. Run this script when:
#   - Setting up a new development environment
#   - After cloning the repository
#   - When you've forgotten which packages you need
#
# Usage:
#   source("dev-setup.R")
#
# =============================================================================

cat("=============================================================================\n")
cat("Installing Development Dependencies for selection.index Package\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Package Development Tools
# -----------------------------------------------------------------------------
cat("1. Installing Package Development Tools...\n")

dev_tools <- c(
  "devtools", # Package development workflow
  "usethis", # Automation for package development tasks
  "roxygen2", # Documentation generation from inline comments
  "pkgdown", # Build package website
  "rcmdcheck" # Run R CMD check from R
)

for (pkg in dev_tools) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   Installing:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("   Already installed:", pkg, "\n")
  }
}

# -----------------------------------------------------------------------------
# 2. Testing Framework
# -----------------------------------------------------------------------------
cat("\n2. Installing Testing Framework...\n")

testing_pkgs <- c(
  "testthat", # Unit testing framework
  "covr" # Code coverage calculation
)

for (pkg in testing_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   Installing:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("   Already installed:", pkg, "\n")
  }
}

# -----------------------------------------------------------------------------
# 3. Documentation & Vignettes
# -----------------------------------------------------------------------------
cat("\n3. Installing Documentation Tools...\n")

doc_pkgs <- c(
  "knitr", # Dynamic report generation
  "rmarkdown", # R Markdown document conversion
  "markdown" # Markdown rendering for R
)

for (pkg in doc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   Installing:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("   Already installed:", pkg, "\n")
  }
}

# -----------------------------------------------------------------------------
# 4. Runtime Dependencies (from DESCRIPTION)
# -----------------------------------------------------------------------------
cat("\n4. Verifying Runtime Dependencies...\n")

runtime_pkgs <- c(
  "stats", # Statistical functions (qt, pf)
  "utils" # Utility functions (combn)
)

for (pkg in runtime_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   Installing:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("   Already installed:", pkg, "\n")
  }
}

# -----------------------------------------------------------------------------
# 5. Optional: Code Quality Tools
# -----------------------------------------------------------------------------
cat("\n5. Installing Code Quality Tools (Optional)...\n")

quality_pkgs <- c(
  "lintr", # Static code analysis
  "styler", # Code formatting
  "goodpractice" # Best practices checker for R packages
)

for (pkg in quality_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   Installing:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat("   Already installed:", pkg, "\n")
  }
}

# -----------------------------------------------------------------------------
# 6. Install Package Dependencies from DESCRIPTION
# -----------------------------------------------------------------------------
cat("\n6. Installing Package Dependencies from DESCRIPTION...\n")

if (requireNamespace("remotes", quietly = TRUE)) {
  remotes::install_deps(dependencies = TRUE)
} else {
  cat("   Installing remotes first...\n")
  install.packages("remotes")
  remotes::install_deps(dependencies = TRUE)
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=============================================================================\n")
cat("Development Environment Setup Complete!\n")
cat("=============================================================================\n\n")

cat("Next steps:\n")
cat("  1. Load the development package: devtools::load_all()\n")
cat("  2. Run tests: devtools::test()\n")
cat("  3. Check package: devtools::check()\n")
cat("  4. Build documentation: devtools::document()\n")
cat("  5. Build website: pkgdown::build_site()\n\n")

cat("Useful commands:\n")
cat("  - Install package locally: devtools::install()\n")
cat("  - Run specific test file: testthat::test_file('tests/testthat/test-*.R')\n")
cat("  - Check code coverage: covr::package_coverage()\n")
cat("  - Run code linting: lintr::lint_package()\n")
cat("  - Check best practices: goodpractice::gp()\n\n")

cat("=============================================================================\n")
