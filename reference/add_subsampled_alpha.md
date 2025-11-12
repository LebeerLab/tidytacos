# Add alpha diversity measures using subsampling

`add_subsampled_alpha()` adds selected alpha diversity measures to the
sample table of a tidytacos object using an itterative subsampling
process.

## Usage

``` r
add_subsampled_alpha(
  ta,
  min_lib_size = NULL,
  method = "shannon",
  itterations = 100
)
```

## Arguments

- ta:

  a tidytacos object.

- min_lib_size:

  the minimum lib size samples need to have. Samples with lower lib
  sizes will be discarded and samples with a higher readcount will be
  itteratively subsampled to this readcount to allow for a fair
  comparison across read_depths.

- method:

  The diversity measure to use, see
  [`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)
  for further information on these.

- itterations:

  the amount of itterations for subsampling. Please report this number
  in your research.

## Value

A tidytacos object with the selected alpha diversity measure added.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

Other diversity-metrics:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md)
