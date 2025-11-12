# Return an interactive ordination plot of the samples

Creates an interactive ordination plot of the beta diversity of the
samples in the tidytacos object. This can be used to gauge the
similarity between samples.

## Usage

``` r
tacoplot_ord_ly(
  ta,
  x = NULL,
  samplenames = sample_id,
  ord = "pcoa",
  dims = 2,
  distance = "bray",
  stat.method = NULL,
  palette = NULL,
  title = NULL,
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- x:

  A string, representing the column name used to color the sample groups
  on.

- samplenames:

  the column in the sample table with the samplenames, defaults to
  sample_id.

- ord:

  the ordination technique to use. Choice from pcoa, tsne and umap.

- dims:

  the amount of dimensions to plot, 2 or 3.

- distance:

  the distance algorithm to use, see
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html).

- stat.method:

  the statistic to print on the figure, choice from mantel and anosim.

- palette:

  A vector of colors, used as the palette for coloring sample

- title:

  a string to display as title of the plot. groups.

- ...:

  Extra arguments to pass to the add_ord function.
