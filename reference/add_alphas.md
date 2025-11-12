# Add alpha diversity measures

[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md)
adds selected alpha diversity measures to the sample table of a
tidytacos object.

## Usage

``` r
add_alphas(ta, methods = "all", ...)
```

## Arguments

- ta:

  A tidytacos object.

- methods:

  A character vector of the diversity measure to use, see
  [`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md)
  for examples. Optionally use 'all' to add all diversity measures.

- ...:

  Arguments passed on to
  [`add_alpha`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md)

  `method`

  :   The diversity measure to use, see
      [`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)
      for further information on these.

  `keep_empty_samples`

  :   Whether to keep empty samples or not.

  `subsample`

  :   Wether to use subsampling to be able to compare samples of varying
      sequencing depths. [Schloss et al.,
      2023](https://journals.asm.org/doi/10.1128/msphere.00355-23)

## Value

A tidytacos object with the selected alpha diversity measures added.

## Details

This function can add multiple different alpha diversity measures to the
sample table, specified by the methods argument.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

Other diversity-metrics:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md)

## Examples

``` r
urt_all_alphas <- urt %>% add_alphas()
#> Warning: Removed 3 empty samples.
```
