# Return an interactive bar plot of the samples

Plots an interactive stacked bar plot of the samples in the tidytacos
object to inspect the taxonomic profile.

## Usage

``` r
tacoplot_stack_ly(ta, n = 12, x = sample_clustered, order_by = NULL)
```

## Arguments

- ta:

  A tidytacos object.

- n:

  An integer, representing the amount of colors used to depict different
  taxa.

- x:

  A string, representing the column name used to label the x-axis

- order_by:

  an optional column name to order the samples by. For examples
  order_by=sample would order the x-axis by the sample names instead of
  by similar profiles.
