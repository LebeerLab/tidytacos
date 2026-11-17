# Return a silhouette score plot to visualize how well PAM clusters are defined, for each number of clusters k in the k range.

Return a silhouette score plot to visualize how well PAM clusters are
defined, for each number of clusters k in the k range.

## Usage

``` r
tacoplot_silhouette_width(ta, k_range = 2:10)
```

## Arguments

- ta:

  A tidytacos object.

- k_range:

  range of number of clusters to divide the data in.

## Examples

``` r
urt %>% tacoplot_silhouette_width(2:6)
#> Warning: No phylogenetic tree found, generating tree using msa of the 'sequence' column.
#> Determining distance matrix based on shared 8-mers:
#> ================================================================================
#> 
#> Time difference of 15.94 secs
#> 
#> Clustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.37 secs
#> 
#> Aligning Sequences:
#> ================================================================================
#> 
#> Time difference of 7.72 secs
#> 
#> Iteration 1 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 1.67 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.42 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 4.37 secs
#> 
#> Iteration 2 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 1.61 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.45 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.96 secs
#> 
#> Refining the alignment:
#> ================================================================================
#> 
#> Time difference of 0.11 secs
#> 
#> 
#> Attaching package: ‘phyloseq’
#> The following objects are masked from ‘package:tidytacos’:
#> 
#>     filter_taxa, rank_names
```
