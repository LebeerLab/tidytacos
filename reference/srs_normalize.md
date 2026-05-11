# Normalize counts using Scaling with Ranked Subsampling (SRS) `srs_normalize()` uses the SRS method to equalize sampling depth to a chosen value.

Normalize counts using Scaling with Ranked Subsampling (SRS)
`srs_normalize()` uses the SRS method to equalize sampling depth to a
chosen value.

## Usage

``` r
srs_normalize(ta, target_reads)
```

## Arguments

- ta:

  A tidytacos object.

- target_reads:

  the sequencing depth to which samples are subsampled. Samples below
  this threshold are discarded.

## Value

A tidytacos object where the reads are subsampled to the chosen value
using SRS.

## See also

Other count-modifiers:
[`add_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_absolute_abundance.md),
[`add_density()`](https://lebeerlab.github.io/tidytacos/reference/add_density.md),
[`add_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_rel_abundance.md)
