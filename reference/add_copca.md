# Add compositional principal components to the sample table

`add_copca()` performs a principal components analysis and adds the
first two principal components to the sample table under column names
"pca_1" and "pca_2".

## Usage

``` r
add_copca(ta)
```

## Arguments

- ta:

  A tidytacos object.

## Value

A tidytacos object with the first two PCA dimensions added to the sample
table.

## Details

Note that this function uses only the 50 most prevalant taxa unless
[`add_logratio()`](https://lebeerlab.github.io/tidytacos/reference/add_logratio.md)
was executed with another value for 'max_taxa'.
