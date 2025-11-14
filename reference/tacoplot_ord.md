# Return an ordination plot of the samples

Creates an ordination plot of the beta diversity of the samples in the
tidytacos object. This can be used to gauge the similarity between
samples.

## Usage

``` r
tacoplot_ord(
  ta,
  x = NULL,
  palette = NULL,
  ord = "pcoa",
  distance = "bray",
  stat.method = NULL,
  title = NULL,
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- x:

  The column name used to color the sample groups on.

- palette:

  A vector of colors, used as the palette for coloring sample groups.

- ord:

  the ordination technique to use. Choice from pcoa, tsne and umap.

- distance:

  the distance algorithm to use, see
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html).

- stat.method:

  the statistic to print on the figure, choice from mantel and anosim.

- title:

  a string to display as title of the plot.

- ...:

  Extra arguments to pass to the add_ord function.

## Examples

``` r
tacoplot_ord(urt, x = location)
#> Warning: Removed 3 empty samples.


# set plate to character, to avoid it being treated as a continuous variable
urt <- urt %>% mutate_samples(plate = as.character(plate))
tacoplot_ord(urt, x = plate, ord = "umap", distance = "aitchison", stat.method = "permanova")
#> Warning: Removed 3 empty samples.
#> Warning: Using pseudocount of 1

```
