# Performs SparCC network analysis on a tidytacos object and then performs Markov Clustering on the network to annotate taxa of the largest clusters in the tidytacos object.

Performs SparCC network analysis on a tidytacos object and then performs
Markov Clustering on the network to annotate taxa of the largest
clusters in the tidytacos object.

## Usage

``` r
cluster_taxa(
  ta,
  min_occurrence = 0.05,
  network_thresh = 0.1,
  min_n_cluster = 3,
  taxon_name = taxon,
  sample_name = sample
)
```

## Arguments

- ta:

  a tidytacos object.

- min_occurrence:

  Percentage of samples the taxon needs to be present in for it to be
  considered in the analysis.

- network_thresh:

  absolute value of correlations below this threshold are filtered out.

- min_n_cluster:

  minimum number of taxa per cluster, smaller clusters are filtered out.

- taxon_name:

  unique name to use for the taxa, by default taxon_id is used.

- sample_name:

  unique name to use for the samples, by default sample_id is used.
