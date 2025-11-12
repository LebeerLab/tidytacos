# Get prevalences of taxa in general or per condition

Returns a tidy table of prevalences: taxon presence counts in samples,
overall or per condition.

## Usage

``` r
prevalences(ta, condition = NULL, pres_abs = F)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  A string denoting a categorical variable in the sample table.

- pres_abs:

  Whether to resort to presence/absence screening.

## Details

Condition should be a categorical variable present in the samples table.
Supply condition as a string.
