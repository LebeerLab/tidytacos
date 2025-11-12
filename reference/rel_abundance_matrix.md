# Return a relative abundance matrix

This function returns a relative abundance matrix; the rows are samples
and the column are taxa.

## Usage

``` r
rel_abundance_matrix(ta, sample_name = sample_id, taxon_name = taxon_id)
```

## Arguments

- ta:

  A tidytacos object.

- sample_name:

  The name of the variable in the sample table to use as row names
  (unquoted).

- taxon_name:

  The name of the variable in the taxon table to use as column names
  (unquoted).

## Value

A matrix with abundance values.
