# Convert matrix with counts to tidy data frame

`counts_tidy()` returns a tidy data frame given a numerical counts
matrix.

## Usage

``` r
counts_tidy(counts_matrix, taxa_are_columns = TRUE, value = "counts")
```

## Arguments

- counts_matrix:

  The count matrix that will be converted.

- taxa_are_columns:

  A logical scalar. Are the taxa defined in columns? Default is TRUE.

- value:

  Name of resulting colum containing the count data. Default is
  "counts".

## Details

This function will convert a numerical counts matrix into a tidy data
frame. To convert a tidy data frame into a numerical counts matrix use
´[`counts_matrix()`](https://lebeerlab.github.io/tidytacos/reference/counts_matrix.md).
