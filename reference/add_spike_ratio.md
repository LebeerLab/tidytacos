# Add spike ratio

`add_spike_ratio()` calculates the ratio of non-spike to spike reads for
each sample and adds this to the sample table under the name
"spike_ratio".

## Usage

``` r
add_spike_ratio(ta, spike_taxon)
```

## Arguments

- ta:

  A tidytacos object.

- spike_taxon:

  The taxon_id of the spike.

## Value

A tidytacos object with the spike ratio added to the sample table.

## Details

This function is useful if a DNA spike was added prior to sequencing and
is based on the method described by [Smets et al.,
2016](https://doi.org/10.1016/j.soilbio.2016.02.003).

Without calculating absolute abundances, the spike ratio allows to
compare absolute abundances between sample. For example, if the spike
ration of one sample is twice that of another, then the absolute number
of sequenced strands at the time of spiking in the one sample is twice
that of the other sample.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

## Examples

``` r
# Initiate counts matrix
x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x,
  taxa_are_columns = FALSE
)

# Add total abundance
data <- data %>%
  add_spike_ratio(spike_taxon = "t1")
```
