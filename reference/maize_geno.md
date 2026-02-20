# Synthetic Maize Genomic Data

A synthetic dataset containing Single Nucleotide Polymorphism (SNP)
marker data for the 100 maize genotypes found in maize_pheno. Designed
for testing genomic selection indices (GESIM/RGESIM) and relationship
matrices.

## Usage

``` r
data(maize_geno)
```

## Format

A numeric matrix with 100 rows (genotypes) and 500 columns (SNP
markers):

- Rows:

  100 genotypes, matching the Genotype column in maize_pheno.

- Columns:

  500 simulated SNP markers, coded as 0, 1, or 2.

## Source

Simulated for the selection.index package to provide reproducible
examples.

## Examples

``` r
data(maize_geno)
dim(maize_geno)
#> [1] 100 500
```
