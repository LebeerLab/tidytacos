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

- ...:

  Arguments passed on to
  [`base::sample`](https://rdrr.io/r/base/sample.html)

  `x`

  :   either a vector of one or more elements from which to choose, or a
      positive integer. See ‘Details.’

  `size`

  :   a non-negative integer giving the number of items to choose.

  `prob`

  :   a vector of probability weights for obtaining the elements of the
      vector being sampled.

## Value

A tidytacos object.
