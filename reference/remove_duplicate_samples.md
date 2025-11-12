# Removes duplicate samples from the tidytacos object

Will remove rows with the exact same metadata as another row but a
different sample_id.

## Usage

``` r
remove_duplicate_samples(ta)
```

## Arguments

- ta:

  a tidytacos object

## Value

the tidytacos object minus the duplicate samples
