# Clusters samples into n clusters

`cluster_samples()` clusters the samples into n clusters and adds these
clusters to a new variable "cluster" in the sample table.

## Usage

``` r
cluster_samples(ta, n_clusters)
```

## Arguments

- ta:

  A tidytacos object.

- n_clusters:

  The number of desired clusters.

## Value

A tidytacos object with a cluster column in the samples table.

## Details

This function calculates the Bray-Curtis distance between samples
followed by hierarchical average linkage clustering of samples. The user
provides a number of desired clusters which will be used to assign the
samples to. A new variable named "cluster" will be added to the samples
tibble of a tidytacos object defining to what cluster a sample belongs.

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
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md)

## Examples

``` r
# Initiate count matrix
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
  cluster_samples(n_clusters = 2)
```
