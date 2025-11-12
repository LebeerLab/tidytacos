# Calculates eigentaxa values based on SparCC - MCL generated clusters per sample. It is advised to run [`cluster_taxa()`](https://lebeerlab.github.io/tidytacos/reference/cluster_taxa.md) on the tidytacos object before running this function to add the clusters if you want to stray from any default parameters.

Calculates eigentaxa values based on SparCC - MCL generated clusters per
sample. It is advised to run
[`cluster_taxa()`](https://lebeerlab.github.io/tidytacos/reference/cluster_taxa.md)
on the tidytacos object before running this function to add the clusters
if you want to stray from any default parameters.

## Usage

``` r
add_eigentaxa(ta, taxon_name = taxon, sample_name = sample)
```

## Arguments

- ta:

  a tidytacos object.

- taxon_name:

  Column name of the taxon identifier, by default taxon.

- sample_name:

  Column name of the sample identifier, by default sample.
