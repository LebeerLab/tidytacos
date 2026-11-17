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
#> Time difference of 12.78 secs
#> 
#> Clustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.28 secs
#> 
#> Aligning Sequences:
#> ================================================================================
#> 
#> Time difference of 5.59 secs
#> 
#> Iteration 1 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 1.29 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.34 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 4.23 secs
#> 
#> Iteration 2 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 1.24 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.33 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.65 secs
#> 
#> Refining the alignment:
#> ================================================================================
#> 
#> Time difference of 0.08 secs
#> 
#> 
#> Attaching package: ‘phyloseq’
#> The following objects are masked from ‘package:tidytacos’:
#> 
#>     filter_taxa, rank_names
```
