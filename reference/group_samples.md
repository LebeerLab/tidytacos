# Group the samples

Group the samples

## Usage

``` r
group_samples(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Grouping criteria for the samples table.

## Value

A (named) list of tidytacos object.

## Examples

``` r
urt_by_loc <- urt %>% group_samples(location)

# apply a function to each separate taco, eg. tacosum
urt_by_loc@tacos %>% lapply(tacosum)
#> $N
#> n_samples    n_taxa   n_reads 
#>        86      1197   1852959 
#> 
#> $NF
#> n_samples    n_taxa   n_reads 
#>       131      1244   2020519 
#> 

# subset urt to keep only nasopharynx samples
urt_nf <- urt_by_loc@tacos$NF
```
