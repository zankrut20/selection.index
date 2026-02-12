# Extract Vector Elements

Extracts specific rows from a column of a matrix. Used for extracting
weight vectors for trait combinations.

## Usage

``` r
cpp_extract_vector(mat, row_indices, col_index)
```

## Arguments

- mat:

  Numeric matrix

- row_indices:

  Integer vector of row indices (1-based)

- col_index:

  Column index (0-based)

## Value

Vector of extracted elements
