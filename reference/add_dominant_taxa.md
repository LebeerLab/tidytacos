# Adds the most dominant taxa in the sample to the sample table

`add_dominant_taxa()` adds the most dominant taxa in the sample above a
given relative abundance to the sample table. Samples that do not have a
dominant taxon will have NA in the dominant_taxon column.

## Usage

``` r
add_dominant_taxa(ta, threshold_dominance = 0.5, taxon_name = taxon_id)
```

## Arguments

- ta:

  A tidytacos object.

- threshold_dominance:

  The relative abundance threshold for a taxon to be considered
  dominant.

- taxon_name:

  The column name of the taxa table that defines the taxon name.

## Value

A tidytacos object with the dominant taxa added to the sample table and
the Berger-Parker index of the most dominant taxon. If the B-P index is
lower than the threshold NA is returned.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)
