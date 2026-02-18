# Phenomic-Genomic Variance-Covariance Matrix (Φ)

Computes the combined phenomic-genomic variance-covariance matrix (Φ or
P_L), which is the block matrix representing the joint distribution of
phenotypes and GEBVs.

Structure: Φ = \[\[P, P_yγ\], \[P_yγ', Γ\]\]

where: - P = Var(y) = phenotypic variance-covariance - Γ = Var(γ) =
genomic variance-covariance - P_yγ = Cov(y, γ) = covariance between
phenotypes and GEBVs

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

  Optional. Phenotypic variance-covariance matrix (n_traits x n_traits).
  If NULL, computed from phen_mat.

- Gamma:

  Optional. Genomic variance-covariance matrix (n_traits x n_traits). If
  NULL, computed from gebv_mat.

- P_yg:

  Optional. Covariance between phenotypes and GEBVs (n_traits x
  n_traits). If NULL, computed from phen_mat and gebv_mat.

- method:

  Method for covariance estimation: "pearson" (default), "kendall",
  "spearman"

- use:

  How to handle missing values: "everything", "all.obs", "complete.obs"
  (default), "na.or.complete", "pairwise.complete.obs"

## Value

Combined phenomic-genomic variance-covariance matrix (2t x 2t) where t
is the number of traits

Combined phenomic-genomic variance-covariance matrix (2t x 2t) where t
is the number of traits

## Details

The phenomic-genomic matrix combines phenotypic and genomic information
in a single covariance structure. It is used in: - CLGSI (Combined
Linear Genomic Selection Index) - GESIM (Genomic Eigen Selection Index
Method) - PPG-GESIM (Predetermined Proportional Gains GESIM)

When GEBVs are highly accurate predictors of true breeding values, P_yγ
≈ Γ, but empirical estimation is more robust for real data.

The phenomic-genomic matrix combines phenotypic and genomic information
in a single covariance structure. It is used in: - CLGSI (Combined
Linear Genomic Selection Index) - GESIM (Genomic Eigen Selection Index
Method) - PPG-GESIM (Predetermined Proportional Gains GESIM)

When GEBVs are highly accurate predictors of true breeding values, P_yγ
≈ Γ, but empirical estimation is more robust for real data.

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapters 4 &
8.

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapters 4 &
8.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate example data
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate phenotypes and GEBVs
set.seed(123)
n_genotypes <- 100
n_traits <- 7
phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                   nrow = n_genotypes, ncol = n_traits)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
                   nrow = n_genotypes, ncol = n_traits)

# Compute phenomic-genomic covariance
Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)
print(dim(Phi))  # Should be 14 x 14 (2 * 7 traits)
} # }
if (FALSE) { # \dontrun{
# Generate example data
gmat <- gen_varcov(seldata[,3:9], seldata[,2], seldata[,1])
pmat <- phen_varcov(seldata[,3:9], seldata[,2], seldata[,1])

# Simulate phenotypes and GEBVs
set.seed(123)
n_genotypes <- 100
n_traits <- 7
phen_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 15, sd = 3),
                   nrow = n_genotypes, ncol = n_traits)
gebv_mat <- matrix(rnorm(n_genotypes * n_traits, mean = 10, sd = 2),
                   nrow = n_genotypes, ncol = n_traits)

# Compute phenomic-genomic covariance
Phi <- phenomic_genomic_varcov(phen_mat, gebv_mat)
print(dim(Phi))  # Should be 14 x 14 (2 * 7 traits)
} # }
```
