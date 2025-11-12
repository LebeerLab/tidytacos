# Add clustering-based sample order

`add_sample_clustered()` adds a new variable defining a sample order
based on a hierarchical clustering of the samples.

## Usage

``` r
add_sample_clustered(ta)
```

## Arguments

- ta:

  A tidytacos object.

## Value

A tidytacos object with a new variable `sample_clustered` added to the
sample table.

## Details

This function calculates the Bray-Curtis distances between samples
followed by hierarchical average linkage clustering of samples. It will
then add a new factor variable "sample_clustered" to the sample tibble
of a tidytacos object. This function is useful if one wants to plot
similar samples together.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

## Examples

``` r
urtc <- urt %>% add_sample_clustered()
urtc$samples %>% dplyr::select(sample_id, sample_clustered)
#> # A tibble: 217 × 2
#>    sample_id sample_clustered
#>    <chr>     <fct>           
#>  1 s1        s1              
#>  2 s2        s2              
#>  3 s3        s3              
#>  4 s4        s4              
#>  5 s5        s5              
#>  6 s6        s6              
#>  7 s7        s7              
#>  8 s8        s8              
#>  9 s9        s9              
#> 10 s10       s10             
#> # ℹ 207 more rows
```
