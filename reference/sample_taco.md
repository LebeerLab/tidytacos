# Take a sample of a tidytacos object

`sample_taco()` takes a random subset of samples from the original
tidytacos object.

## Usage

``` r
sample_taco(
  tt,
  n,
  group_col,
  taxon_identifier = sequence,
  replace = FALSE,
  ...
)
```

## Arguments

- tt:

  A tidytacos object.

- n:

  The amount of samples that need to be returned.

- group_col:

  An optional name of a field in the sample table which the samples need
  to be evenly distributed over.

- taxon_identifier:

  The column that uniquely identifies a taxa

- replace:

  Replace selected samples so they can be picked again in sampling.

## Value

A tidytacos object.
