# Genomic Variance-Covariance Matrix (Γ)

Computes genomic variance-covariance matrix (Γ or Gamma) from a matrix
of Genomic Estimated Breeding Values (GEBVs).

γ (gamma) represents GEBV vectors obtained from genomic prediction
models (e.g., GBLUP, rrBLUP, Genomic BLUP). This function computes
Var(γ) = Γ.

## Usage

``` r
genomic_varcov(gebv_mat, method = "pearson", use = "complete.obs")
```

## Arguments

- gebv_mat:

  Matrix of GEBVs (n_genotypes x n_traits)

- method:

  Character string specifying correlation method: "pearson" (default),
  "kendall", or "spearman"

- use:

  Character string specifying how to handle missing values: "everything"
  (default), "complete.obs", "pairwise.complete.obs", etc. See
  [`cov`](https://rdrr.io/r/stats/cor.html) for details.

## Value

Symmetric genomic variance-covariance matrix (n_traits x n_traits)

## Details

The genomic variance-covariance matrix Γ captures genetic variation as
predicted by molecular markers. It is computed as:

where γ_i is the GEBV vector for genotype i and μ_γ is the mean GEBV
vector.

\*\*Missing Value Handling:\*\* - "complete.obs": Uses only complete
observations (recommended) - "pairwise.complete.obs": Uses
pairwise-complete observations (may not be PSD) - "everything": Fails if
any NA present

When using pairwise deletion, the resulting matrix may not be positive
semi-definite (PSD), which can cause numerical issues in selection
indices.

\*\*Applications:\*\* In selection index theory: - Used in LGSI (Linear
Genomic Selection Index) - Component of Φ (phenomic-genomic
covariance) - Component of A (genetic-genomic covariance)

## References

Cerón-Rojas, J. J., & Crossa, J. (2018). Linear Selection Indices in
Modern Plant Breeding. Springer International Publishing. Chapters 4 &
8.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate GEBVs
set.seed(123)
n_genotypes <- 100
n_traits <- 5
gebv_mat <- matrix(rnorm(n_genotypes * n_traits),
  nrow = n_genotypes, ncol = n_traits
)
colnames(gebv_mat) <- paste0("Trait", 1:n_traits)

# Compute genomic variance-covariance
Gamma <- genomic_varcov(gebv_mat)
print(Gamma)
} # }
```
