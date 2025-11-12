# Return a list of taxon IDs per condition

This function returns a named list of taxon_ids per distinct value of a
categorical column of the samples table.

## Usage

``` r
taxonlist_per_condition(ta, condition, read_treshold = 0)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  The name of a variable in the sample table that contains a categorical
  value.

- read_treshold:

  The minimum read count to consider a taxon.

## Value

A list of taxon_id vectors.
