# Return a counts matrix

This function returns a matrix with taxon counts; the rows are samples
and the columns are taxa.

## Usage

``` r
counts_matrix(
  ta,
  sample_name = sample_id,
  taxon_name = taxon_id,
  value = count,
  keep_empty_samples = FALSE
)
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

- value:

  The name of the variable in the counts table to use as count

- keep_empty_samples:

  Should empty samples be included in the matrix?

## Value

A matrix with count values.
