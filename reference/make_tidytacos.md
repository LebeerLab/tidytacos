# Create a tidytacos object from three tidy tables

Create a tidytacos object from three tidy tables

## Usage

``` r
make_tidytacos(
  samples,
  taxa,
  counts,
  sample_name = sample,
  taxon_name = sequence
)
```

## Arguments

- samples:

  A tidy table containing sample information.

- taxa:

  A tidy table containing taxon information.

- counts:

  A tidy table, where each row represents the counts of a taxon in a
  sample.

- sample_name:

  The column in the sample table that contains a unique identifier for
  each sample.

- taxon_name:

  The column in the taxon table that contains a unique identifier for
  each taxon.
