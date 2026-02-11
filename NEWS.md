# selection.index (development version)

## Phase 1: Enhanced Desired Gains Index (Chapter 4 Implementation)

### New Features in `dg_lpsi()`
* **Implied Economic Weights**: Calculate and return the implied economic weights (ŵ = G^-1 P b) that would be needed in a Smith-Hazel index to achieve the desired gains (Pesek & Baker, 1969, Section 1.4)
  - Added `return_implied_weights` parameter (default: TRUE)
  - Returns both raw and normalized (max = 1) implied weights
  - Weights included in summary data frame

* **Feasibility Checking**: Automatic validation of desired gains against theoretical maximums
  - Added `check_feasibility` parameter (default: TRUE)
  - Added `selection_intensity` parameter to customize feasibility thresholds (default: 2.063 for i = 10% selection)
  - Warns when desired gains exceed 80% of theoretical maximum (i × √G_ii)
  - Returns detailed feasibility metrics including genetic SDs, max possible gains, and feasibility ratios

* **Enhanced Output Structure**: `dg_lpsi()` now returns:
  - `implied_weights`: Vector of implied economic weights
  - `implied_weights_normalized`: Normalized weights (maximum absolute value = 1)
  - `feasibility`: Data frame with per-trait feasibility analysis
  - `gain_errors`: Per-trait precision metrics (desired - achieved)
  - `selection_intensity`: Selection intensity used in calculations

* **S3 Methods**: New print and summary methods for improved output
  - `print.dg_lpsi()`: Formatted display with gains comparison, implied weights, and feasibility analysis
  - `summary.dg_lpsi()`: Additional aggregate statistics including mean gains, feasibility summary, and interpretation guidance
  - Visual indicators (✓, ⚠, ✗) for at-a-glance assessment

### Testing & Documentation
* Added 24 comprehensive unit tests for new functionality
* Created demonstration script: `inst/examples/phase1_demo.R`
* All 97 package tests pass with backward compatibility maintained

### Technical Improvements
* Uses MASS::ginv() for robust matrix inversion in implied weights calculation
* Proper handling of missing trait names in print methods
* Enhanced error messages with mathematical context

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
