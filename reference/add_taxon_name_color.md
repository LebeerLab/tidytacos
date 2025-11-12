# Add taxon color for visualization.

[`add_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_rel_abundance.md)
determines the most abundant taxa and assigns them a color for
consistent color codes of each taxon in visualizations. A rank can be
supplied to aggregate colors higher than the current rank.

## Usage

``` r
add_taxon_name_color(
  ta,
  method = "mean_rel_abundance",
  n = 12,
  samples = NULL,
  taxa = NULL,
  rank = NULL,
  threshold_dominance = NULL
)
```

## Arguments

- ta:

  A tidytacos object.

- method:

  The method by which to arrange the taxon names. Currently only
  mean_rel_abundance or dominance.

- n:

  An integer denoting the amount of most abundant taxa to display.
  Capacity at 12.

- samples:

  An optional vector of sample_id's of interest.

- taxa:

  An optional vector of taxon_id's of interest.

- rank:

  An optional rank to aggregate taxa on.

- threshold_dominance:

  An optional threshold for the dominance method.

## Value

A tidytacos object.

## See also

Other taxa-modifiers:
[`add_mean_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_mean_rel_abundance.md),
[`add_prevalence()`](https://lebeerlab.github.io/tidytacos/reference/add_prevalence.md),
[`add_taxon_name()`](https://lebeerlab.github.io/tidytacos/reference/add_taxon_name.md)

## Examples

``` r
# display the 5 most abundant taxa at genus lvl
urt %>% add_taxon_name_color(n=5, rank='genus') %>% tacoplot_stack()
```
