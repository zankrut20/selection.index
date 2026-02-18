# Genomic Variance-Covariance Matrix (Γ)

Computes the genomic variance-covariance matrix (Γ or G_s) from Genomic
Estimated Breeding Values (GEBVs) or molecular marker scores.

This represents: Γ = Var(γ)

where γ is the vector of GEBVs across traits.

Computes the genomic variance-covariance matrix (Γ or G_s) from Genomic
Estimated Breeding Values (GEBVs) or molecular marker scores.

This represents: Γ = Var(γ)

where γ is the vector of GEBVs across traits.

## Usage

``` r
genomic_varcov(gebv_mat, method = "pearson", use = "complete.obs")

genomic_varcov(gebv_mat, method = "pearson", use = "complete.obs")
```

## Arguments

- gebv_mat:

  Matrix of GEBVs (n_genotypes x n_traits)

- method:

  Method for covariance estimation: "pearson" (default), "kendall",
  "spearman"

- use:

  How to handle missing values: "everything", "all.obs", "complete.obs"
  (default), "na.or.complete", "pairwise.complete.obs"

## Value

Symmetric genomic variance-covariance matrix (n_traits x n_traits)

Symmetric genomic variance-covariance matrix (n_traits x n_traits)

## Details

The genomic variance-covariance matrix represents the
variance-covariance structure of GEBVs predicted from molecular markers.
It captures the genetic relationships among traits as explained by the
markers.

In selection index theory: - Used in LGSI (Linear Genomic Selection
Index) - Component of Φ (phenomic-genomic covariance) - Component of A
(genetic-genomic covariance)

The genomic variance-covariance matrix represents the
variance-covariance structure of GEBVs predicted from molecular markers.
It captures the genetic relationships among traits as explained by the
markers.

In selection index theory: - Used in LGSI (Linear Genomic Selection
Index) - Component of Φ (phenomic-genomic covariance) - Component of A
(genetic-genomic covariance)

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
# Simulate GEBVs
set.seed(123)
n_genotypes <- 100
n_traits <- 5
gebv_mat <- matrix(rnorm(n_genotypes * n_traits), 
                   nrow = n_genotypes, ncol = n_traits)
colnames(gebv_mat) <- paste0("Trait", 1:n_traits)

# Compute genomic variance-covariance
Gamma <- genomic_varcov(gebv_mat)
print(Gamma)
} # }
if (FALSE) { # \dontrun{
# Simulate GEBVs
set.seed(123)
n_genotypes <- 100
n_traits <- 5
gebv_mat <- matrix(rnorm(n_genotypes * n_traits), 
                   nrow = n_genotypes, ncol = n_traits)
colnames(gebv_mat) <- paste0("Trait", 1:n_traits)

# Compute genomic variance-covariance
Gamma <- genomic_varcov(gebv_mat)
print(Gamma)
} # }
```
