# Linear Molecular and Genomic Eigen Selection Index Methods (Chapter 8)

Implements the Linear Molecular and Genomic Eigen Selection Index
methods from Chapter 8. These methods extend eigen-based selection to
genomic/molecular data, maximizing the accuracy squared (rho_HI^2)
through eigenvalue problems without requiring pre-specified economic
weights.

Methods included: - MESIM : Molecular Eigen Selection Index Method
(Section 8.1) - GESIM : Linear Genomic Eigen Selection Index Method
(Section 8.2) - GW-ESIM : Genome-Wide Linear Eigen Selection Index
Method (Section 8.3) - RGESIM : Restricted Linear Genomic Eigen
Selection Index Method (Section 8.4) - PPG-GESIM: Predetermined
Proportional Gain Genomic Eigen Selection Index (Section 8.5)

All implementations use C++ primitives (math_primitives.cpp) for
quadratic forms and symmetric solves, while eigendecompositions use R's
eigen() for correctness and compatibility with the existing package
architecture.

## Mathematical Foundation

Like the phenotypic ESIM (Chapter 7), these genomic eigen methods
maximize the squared accuracy between the index and net genetic merit,
but incorporate molecular markers, GEBVs, or genome-wide marker scores.

The general approach solves a generalized eigenproblem to find optimal
index coefficients that maximize heritability without requiring economic
weights.

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 8.
