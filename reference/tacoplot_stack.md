# Return a bar plot of the samples

Plots a stacked bar plot of the samples in the tidytacos object to
inspect the taxonomic profile.

## Usage

``` r
tacoplot_stack(ta, n = 12, x = sample_clustered, pie = FALSE, order_by = NULL)
```

## Arguments

- ta:

  A tidytacos object.

- n:

  An integer, representing the amount of colors used to depict

- x:

  The name of the column name used to represent samples on the x-axis

- pie:

  A boolean, whether or not to represent the profile in a pie chart.
  Default is FALSE, as pie chart representations can be misleading to
  interpret.

- order_by:

  an optional column name to order the samples by. For examples
  order_by=sample would order the x-axis by the sample names instead of
  by similar profiles.
