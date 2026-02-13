# Constrained Linear Genomic Selection Indices (Chapter 6)

## Overview

This document describes the implementation of constrained linear genomic selection indices from Chapter 6 of Cerón-Rojas & Crossa (2018). All mathematical formulas have been implemented in the new modular file `constrained_genomic_indices.R`.

## Implementation Architecture

### Design Principles

1. **Modular Structure**: All constrained genomic indices in a separate file for easy maintenance
2. **C++ for Math**: All mathematical operations use primitives from `math_primitives.cpp`
3. **R for Logic**: All selection index logic and algorithms implemented in R
4. **Consistent API**: Functions follow the same pattern as existing package functions

### Files Created

- **R Code**: `R/constrained_genomic_indices.R` (1200+ lines)
- **Documentation**: 
  - `man/rlgsi.Rd`
  - `man/ppg_lgsi.Rd`
  - `man/crlgsi.Rd`
  - `man/cppg_lgsi.Rd`
  - `man/constrained_genomic_indices.Rd`
- **Examples**: `inst/examples/chapter6_constrained_genomic_indices.R`
- **Exports**: Updated `NAMESPACE`

---

## Mathematical Formulas → R Implementation

### 6.1 RLGSI (Restricted Linear Genomic Selection Index)

#### LaTeX Formulas (from Chapter 6.tex)

**Objective Function** (Eq. 6.1):
```
f_R(β,v) = w'Cw + β'Γβ - 2w'Γβ + 2v'U'Γβ
```

**Augmented System** (Eq. 6.2):
```
[Γ   ΓU ] [β]   [Γw]
[U'Γ  0 ] [v] = [ 0 ]
```

**Coefficients** (Eq. 6.3):
```
β_RG = K_G w
```

**Selection Response** (Eq. 6.5):
```
R_RG = (k_I / L_G) * sqrt(β_RG' Γ β_RG)
```

**Expected Gains** (Eq. 6.6):
```
E_RG = (k_I / L_G) * (Γ β_RG) / sqrt(β_RG' Γ β_RG)
```

#### R Implementation

```r
rlgsi(Gamma, wmat, wcol = 1, restricted_traits = NULL, U = NULL, 
      k_I = 2.063, L_G = 1, GAY = NULL)
```

**Key Implementation Details**:
- Build augmented matrix using `rbind()`, `cbind()`
- Solve system using `MASS::ginv()` for robustness
- Extract β_RG from solution vector (first n_traits elements)
- Use `cpp_quadratic_form_sym()` for β'Γβ calculation
- Response and gains computed in `.genomic_index_metrics()` helper

---

### 6.2 PPG-LGSI (Predetermined Proportional Gains)

#### LaTeX Formulas

**Alternative Form** (Eq. 6.9):
```
β_PG = β_RG + θ_G * U(U'ΓU)^(-1)d
```

**Proportionality Constant** (Eq. 6.10):
```
θ_G = [d'(U'ΓU)^(-1)U'Γw] / [d'(U'ΓU)^(-1)d]
```

**Selection Response** (Eq. 6.11):
```
R_PG = (k_I / L_G) * sqrt(β_PG' Γ β_PG)
```

**Expected Gains** (Eq. 6.12):
```
E_PG = (k_I / L_G) * (Γ β_PG) / sqrt(β_PG' Γ β_PG)
```

#### R Implementation

```r
ppg_lgsi(Gamma, d, wmat = NULL, wcol = 1, U = NULL,
         k_I = 2.063, L_G = 1, GAY = NULL)
```

**Key Implementation Details**:
- First compute β_RG using internal call to `rlgsi()`
- Compute U'ΓU and invert using `solve()`
- Calculate θ_G using matrix operations
- Compute δ = U(U'ΓU)^(-1)d
- Final coefficients: β_PG = β_RG + θ_G * δ
- Verify proportionality by computing gain ratios

---

### 6.3 CRLGSI (Combined Restricted Linear Genomic Selection Index)

#### LaTeX Formulas

**Coefficients** (Eq. 6.13):
```
β_CR = K_C β_C
```

**Selection Response** (Eq. 6.14):
```
R_CR = (k_I / L_I) * sqrt(β_CR' T_C β_CR)
```

**Expected Gains** (Eq. 6.15):
```
E_CR = (k_I / L_I) * (Ψ_C β_CR) / sqrt(β_CR' T_C β_CR)
```

#### R Implementation

```r
crlgsi(T_C = NULL, Psi_C = NULL, phen_mat = NULL, gebv_mat = NULL,
       pmat = NULL, gmat = NULL, wmat, wcol = 1,
       restricted_traits = NULL, U = NULL, reliability = NULL,
       k_I = 2.063, L_I = 1, GAY = NULL)
```

**Key Implementation Details**:
- Build T_C matrix: [P, P_yg; P_yg', P_g] (2t x 2t)
- Build Ψ_C matrix: [G; C_gebv_g] (2t x t)
- Transform constraint matrix: U_TC = Ψ_C %*% U
- Solve augmented system similar to RLGSI
- Split solution into b_y (phenotype) and b_g (GEBV) coefficients
- Use `.combined_index_metrics()` helper for calculations

---

### 6.4 CPPG-LGSI (Combined Predetermined Proportional Gains)

#### LaTeX Formulas

**Coefficients** (Eq. 6.16):
```
β_CP = β_CR + θ_CP * δ_CP
```

**Proportionality Constant** (Eq. 6.17):
```
θ_CP = [β_C' Φ_C(Φ_C' T_C^(-1) Φ_C)^(-1) d_C] / 
       [d_C' (Φ_C' T_C^(-1) Φ_C)^(-1) d_C]
```

**Selection Response** (Eq. 6.18):
```
R_CP = (k_I / L_I) * sqrt(β_CP' T_C β_CP)
```

**Expected Gains** (Eq. 6.19):
```
E_CP = (k_I / L_I) * (Ψ_C β_CP) / sqrt(β_CP' T_C β_CP)
```

#### R Implementation

```r
cppg_lgsi(T_C = NULL, Psi_C = NULL, d, phen_mat = NULL, gebv_mat = NULL,
          pmat = NULL, gmat = NULL, wmat = NULL, wcol = 1, U = NULL,
          reliability = NULL, k_I = 2.063, L_I = 1, GAY = NULL)
```

**Key Implementation Details**:
- First compute β_CR using internal call to `crlgsi()`
- Compute Φ_C = Ψ_C %*% U
- Invert T_C using `MASS::ginv()`
- Calculate (Φ_C' T_C^(-1) Φ_C)^(-1)
- Compute unrestricted β_C = T_C^(-1) Ψ_C w
- Calculate θ_CP using formula
- Compute δ_CP and final coefficients
- Split into phenotype and GEBV components

---

## C++ Primitives Used

All mathematical operations leverage existing C++ functions from `math_primitives.cpp`:

| Operation | C++ Function | Usage |
|-----------|--------------|-------|
| Solve Ax = b | `cpp_symmetric_solve()` | Solving linear systems |
| x'Ay | `cpp_quadratic_form()` | Computing quadratic forms |
| x'Ax | `cpp_quadratic_form_sym()` | Computing symmetric quadratic forms |
| Extract elements | `cpp_extract_vector()` | Extracting weight vectors |

**No new C++ code was needed** - all existing primitives were sufficient!

---

## Helper Functions

### `.genomic_index_metrics()`

Computes standard metrics for genomic indices (RLGSI, PPG-LGSI):
- Uses Γ (GEBV variance-covariance)
- Computes b'Γb, selection response, expected gains
- Returns accuracy (rHI), overall GA, PRE

### `.combined_index_metrics()`

Computes metrics for combined indices (CRLGSI, CPPG-LGSI):
- Uses T_C (combined variance) and Ψ_C (combined genetic covariance)
- Computes b'T_Cb, selection response, expected gains
- Returns accuracy, GA, PRE

Both helpers use only C++ primitives for calculations.

---

## Usage Examples

### Simple Restriction Example

```r
# Restrict traits 2 and 4 to zero gain
result <- rlgsi(Gamma, weights, restricted_traits = c(2, 4))

# Check constraints (should be ~0)
print(result$E[c(2, 4)])
```

### Proportional Gains Example

```r
# Achieve 3:2:1 ratio for first 3 traits
d <- c(3, 2, 1, 0, 0, 0, 0)
result <- ppg_lgsi(Gamma, d, wmat = weights)

# Verify proportionality
print(result$gain_ratios)  # Should be constant
```

### Combined Data with Restrictions

```r
# Use both phenotypes and GEBVs
result <- crlgsi(
  phen_mat = phenotypes,
  gebv_mat = gebvs,
  pmat = pmat,
  gmat = gmat,
  wmat = weights,
  restricted_traits = c(2, 4),
  reliability = 0.7
)

# Access separate coefficients
print(result$b_y)  # Phenotype coefficients
print(result$b_g)  # GEBV coefficients
```

---

## Testing and Validation

### Constraint Verification

All functions return `constrained_response` to verify constraints are satisfied:
- RLGSI: U'Γβ ≈ 0
- PPG-LGSI: gain_ratios should be constant
- CRLGSI: U'E ≈ 0
- CPPG-LGSI: gain_ratios should be constant

### Numerical Stability

- All functions use `MASS::ginv()` for matrix inversion (handles near-singular matrices)
- Augmented systems with explicit zero blocks for Lagrange multipliers
- Checks for NA/Inf in coefficient vectors

---

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). *Linear Selection Indices in Modern Plant Breeding*. 
Springer International Publishing. Chapter 6: Constrained Linear Genomic Selection Indices.

---

## Package Integration

### Existing Functions

These new functions complement existing package capabilities:

**Unrestricted Indices:**
- `lpsi()` - Linear Phenotypic Selection Index
- `lgsi()` - Linear Genomic Selection Index  
- `clgsi()` - Combined Linear Genomic Selection Index

**Constrained Phenotypic Indices:**
- `rlpsi()` - Restricted LPSI
- `ppg_lpsi()` - Predetermined Proportional Gains LPSI
- `dg_lpsi()` - Desired Gains LPSI

**New Constrained Genomic Indices:**
- `rlgsi()` - Restricted LGSI
- `ppg_lgsi()` - Predetermined Proportional Gains LGSI
- `crlgsi()` - Combined Restricted LGSI
- `cppg_lgsi()` - Combined Predetermined Proportional Gains LGSI

### Consistent API Design

All functions share common parameter names:
- `restricted_traits` - User-friendly way to specify constraints
- `U` - Advanced constraint matrix option
- `k_I` - Selection intensity
- `GAY` - Comparative trait for PRE calculation
- `wmat`, `wcol` - Economic weights

All return lists with:
- `summary` - Data frame with coefficients and metrics
- `b` - Coefficient vector
- `E` - Expected genetic gains per trait
- `R` - Overall selection response

---

## Future Extensions

Possible additions following the same modular pattern:

1. **Print/Summary Methods**: S3 methods for the new index classes
2. **Visualization**: Plot functions for constraint visualization
3. **Cross-validation**: Functions to validate reliability estimates
4. **Multi-trait GEBV**: Support for joint modeling
5. **Economic Weights**: Automated weight estimation from market data

All extensions would follow the same principles:
- C++ for computation
- R for logic
- Modular organization
- Consistent API

---

## Conclusion

The implementation successfully translates all mathematical formulas from Chapter 6 into efficient, well-documented R code. The modular design makes it easy to:

- Understand how formulas map to code
- Maintain and extend functionality
- Test individual components
- Integrate with existing package functions

All code follows package conventions and uses existing C++ primitives, requiring no new C++ development.
