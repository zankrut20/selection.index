# Phenomic-Genomic Variance-Covariance Matrix (Φ)

Computes the combined phenomic-genomic variance-covariance matrix (Φ or
P_L), which is the block matrix representing the joint distribution of
phenotypes and GEBVs.

Structure: Φ = \[\[P, P_yγ\], \[P_yγ', Γ\]\]

where: - P = Var(y) = phenotypic variance-covariance - Γ = Var(γ) =
genomic variance-covariance - P_yγ = Cov(y, γ) = covariance between
phenotypes and GEBVs

## Usage

``` r
phenomic_genomic_varcov(
  phen_mat = NULL,
  gebv_mat = NULL,
  P = NULL,
  Gamma = NULL,
  P_yg = NULL,
  method = "pearson",
  use = "complete.obs"
)
```

## Arguments

- phen_mat:

  Matrix of phenotypes (n_genotypes x n_traits). Optional if P and P_yg
  are provided.

- gebv_mat:

  Matrix of GEBVs (n_genotypes x n_traits). Optional if Gamma and P_yg
  are provided.

- P:

  Phenotypic variance-covariance matrix (n_traits x n_traits). Optional
  if phen_mat is provided.

- Gamma:

  Genomic variance-covariance matrix (n_traits x n_traits). Optional if
  gebv_mat is provided.

- P_yg:

  Covariance between phenotypes and GEBVs (n_traits x n_traits).
  Optional if phen_mat and gebv_mat are provided.

- method:

  Character string specifying correlation method: "pearson" (default),
  "kendall", or "spearman"

- use:

  Character string specifying how to handle missing values:
  "complete.obs" (default), "pairwise.complete.obs", etc.

## Value

Symmetric block matrix Φ (2\*n_traits x 2\*n_traits)

## Details

The phenomic-genomic covariance matrix is used in: - GESIM (Genomic
Eigen Selection Index Method) - Combined phenotypic + genomic selection
indices

The matrix is constructed as: \$\$\Phi = \begin{bmatrix} P &
P\_{y\gamma} \\ P\_{y\gamma}' & \Gamma \end{bmatrix}\$\$

where the off-diagonal blocks are transposes, ensuring symmetry.

You can provide either: 1. Raw data: phen_mat + gebv_mat (matrices
computed internally) 2. Pre-computed matrices: P + Gamma + P_yg

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapter 8.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
n_genotypes <- 100
n_traits <- 7
phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
  nrow = n_genotypes, ncol = n_traits
)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
  nrow = n_genotypes, ncol = n_traits
)

# Compute phenomic-genomic covariance
Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)
print(dim(Phi)) # Should be 14 x 14 (2 * 7 traits)
} # }
```
