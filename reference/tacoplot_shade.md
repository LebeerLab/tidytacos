# Return a bar plot of the samples in the shade of one taxonomic rank; colored by a lower level rank.

Plots a stacked bar plot of the samples in the tidytacos object to
inspect the taxonomic profile on two taxonomic ranks.

## Usage

``` r
tacoplot_shade(
  ta,
  group = "family",
  subgroup = "genus",
  n = 5,
  N = 4,
  order_tax = "NA"
)
```

## Arguments

- ta:

  A tidytacos object.

- group:

  The higher level taxonomic rank name (default: family)

- subgroup:

  The lower level taxonomic rank name (default: genus)

- n:

  Number of groups to plot.

- N:

  Number of subgroups per group; max of 4.

- order_tax:

  an optional individual name of the group taxonomic rank to sort the
  samples on. For examples order_by=Lactobacillaceae would order the
  x-axis by the abundance of Lactobacillaceae.
