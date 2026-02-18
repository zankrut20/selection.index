# Linear Phenotypic Eigen Selection Index Methods (Chapter 7)

Implements the Linear Phenotypic Eigen Selection Index methods from
Chapter 7. These methods resolve index coefficients by maximizing the
accuracy squared (rho_HI^2) through an eigenvalue problem rather than
requiring pre-specified economic weights.

Methods included: - ESIM : Linear Phenotypic Eigen Selection Index
(Section 7.1) - RESIM : Linear Phenotypic Restricted Eigen Selection
Index (Section 7.2) - PPG-ESIM: Predetermined Proportional Gain Eigen
Selection Index (Section 7.3)

All implementations use C++ primitives (math_primitives.cpp) for
quadratic forms and symmetric solves, while eigendecompositions use R's
eigen() for correctness and compatibility with the existing package
architecture.

## Mathematical Foundation

Unlike classical LPSI which requires economic weights w, the ESIM family
resolves the index vector b_E by maximizing the squared accuracy:
\$\$\rho\_{HI}^2 = \frac{b'Cb}{b'Pb}\$\$ leading to the generalized
eigenproblem \\(\mathbf{P}^{-1}\mathbf{C} - \lambda^2 I)b = 0\\.

The largest eigenvalue lambda_E^2 equals the maximum achievable index
heritability, and the corresponding eigenvector b_E contains the optimal
index coefficients.

## References

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 7.

Ceron-Rojas, J. J., Crossa, J., Sahagun-Castellanos, J.,
Castillo-Gonzalez, F., & Santacruz-Varela, A. (2006). A selection index
method based on eigen analysis. Crop Science, 46(4), 1711-1721.
