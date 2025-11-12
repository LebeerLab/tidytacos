# Add total densities of samples

`add_total_density()` adds the total microbial density to the sample
table of a tidytacos object under the column name "total_density".

## Usage

``` r
add_total_density(
  ta,
  spike_taxon,
  spike_added = spike_added,
  material_sampled = material_sampled
)
```

## Arguments

- ta:

  A tidytacos object.

- spike_taxon:

  The taxon id of the spike.

- spike_added:

  The column name of the samples table which indicates how much spike
  was added per sample, e.g. 16S rRNA gene copy numbers added to the DNA
  extraction tube.

- material_sampled:

  The column name indicating the amount of material from which DNA was
  extracted, e.g gram of soil. This parameter encourages researchers to
  consider that absolute abundances are only meaningful if they can be
  translated into densities.

## Value

A tidytacos object with the total densities added to the sample table.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

## Examples

``` r
# Initiate count matrix
x <- matrix(
  c(1500, 1300, 14, 280, 356, 9),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2", "taxon3")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x,
  taxa_are_columns = FALSE
)
data$samples$spike_added <- c(100, 150)
data$samples$material_sampled <- c(1, 5)

# Add total abundance
data <- data %>%
  add_total_density(spike_taxon = "t3")
```
