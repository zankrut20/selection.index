# Generic C++ Math Primitives for Experimental Design Statistics

Generic mathematical operations optimized with C++/Eigen. No
design-specific logic - purely mathematical primitives that can be
orchestrated by R code to implement any experimental design.

This architecture allows: - Easy addition of new experimental designs
(in R only) - C++ speed for heavy computation - Single source of truth
(design_stats.R) - Better maintainability and testability
