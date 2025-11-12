# Perform anosim test

`perform_anosim()` performs the anosim test for statistical difference
between groups of samples. The null hypothesis is that there is no
difference between microbial communities in the groups of samples.

## Usage

``` r
perform_anosim(ta, group, ...)
```

## Arguments

- ta:

  A tidytacos object.

- group:

  A column in the sample table to group the samples on.

- ...:

  Arguments passed on to
  [`vegan::anosim`](https://vegandevs.github.io/vegan/reference/anosim.html)

  `x`

  :   Data matrix or data frame in which rows are samples and columns
      are response variable(s), or a dissimilarity object or a symmetric
      square matrix of dissimilarities.

  `grouping`

  :   Factor for grouping observations.

  `permutations`

  :   a list of control values for the permutations as returned by the
      function [`how`](https://rdrr.io/pkg/permute/man/how.html), or the
      number of permutations required, or a permutation matrix where
      each row gives the permuted indices.

  `distance`

  :   Choice of distance metric that measures the dissimilarity between
      two observations. See
      [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
      for options. This will be used if `x` was not a dissimilarity
      structure or a symmetric square matrix.

  `strata`

  :   An integer vector or factor specifying the strata for permutation.
      If supplied, observations are permuted only within the specified
      strata.

  `parallel`

  :   Number of parallel processes or a predefined socket cluster. With
      `parallel = 1` uses ordinary, non-parallel processing. The
      parallel processing is done with parallel package.

## Value

an object of class "anosim"

## Examples

``` r
perform_anosim(urt, method, dist = "jaccard")
#> Warning: Empty samples found, ignoring them in analysis
#> Warning: Removed 3 empty samples.
#> 
#> Call:
#> vegan::anosim(x = M, grouping = ta$samples %>% pull(!!group),      distance = "jaccard") 
#> Dissimilarity: jaccard 
#> 
#> ANOSIM statistic R: -0.05501 
#>       Significance: 0.827 
#> 
#> Permutation: free
#> Number of permutations: 999
#> 
# no statistical difference based on the method column
# (high significance value and R close to 0).
```
