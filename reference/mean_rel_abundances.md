# Get mean relative abundances of taxa in general or per condition

Returns tidy table with average relatively abundances of taxa, overall
or per condition.

## Usage

``` r
mean_rel_abundances(ta, condition = NULL)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  A string representing a categorical variable to compute the relative
  abundances in every option of the variable

## Details

Condition should be a categorical variable present in the samples table.
Supply condition as a string.
