# Aggregate samples with identical values for all metadata

`aggregate_samples()` merges sample content of samples which have
identical values for all columns in the sample table (except sample_id).

## Usage

``` r
aggregate_samples(ta)
```

## Arguments

- ta:

  A tidytacos object.

## Value

A tidytacos object.
