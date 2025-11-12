# Perform an adonis test

This function executes the
[adonis2](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)
function of the vegan package and returns the result.

## Usage

``` r
perform_adonis(ta, predictors, permutations = 999, ...)
```

## Arguments

- ta:

  A tidytacos object.

- predictors:

  A character vector with predictors to include in the model.

- permutations:

  The number of permutations (more permutations takes longer but gives a
  more accurate p-value).

- ...:

  Arguments passed on to
  [`vegan::adonis2`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)

  `formula`

  :   Model formula. The left-hand side (LHS) of the formula must be
      either a community data matrix or a dissimilarity matrix, e.g.,
      from
      [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
      or [`dist`](https://rdrr.io/r/stats/dist.html). If the LHS is a
      data matrix, function
      [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
      will be used to find the dissimilarities. The right-hand side
      (RHS) of the formula defines the independent variables. These can
      be continuous variables or factors, they can be transformed within
      the formula, and they can have interactions as in a typical
      [`formula`](https://rdrr.io/r/stats/formula.html).

  `data`

  :   the data frame for the independent variables, with rows in the
      same order as the community data matrix or dissimilarity matrix
      named on the LHS of `formula`.

  `method`

  :   the name of any method used in
      [`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
      to calculate pairwise distances if the left hand side of the
      `formula` was a data frame or a matrix.

  `sqrt.dist`

  :   Take square root of dissimilarities. This often euclidifies
      dissimilarities.

  `add`

  :   Add a constant to the non-diagonal dissimilarities such that all
      eigenvalues are non-negative in the underlying Principal
      Co-ordinates Analysis (see
      [`wcmdscale`](https://vegandevs.github.io/vegan/reference/wcmdscale.html)
      for details). Choice `"lingoes"` (or `TRUE`) use the recommended
      method of Legendre & Anderson (1999: “method 1”) and `"cailliez"`
      uses their “method 2”.

  `by`

  :   `by = NULL` will assess the overall significance of all terms
      together, `by = "terms"` will assess significance for each term
      (sequentially from first to last), setting `by = "margin"` will
      assess the marginal effects of the terms (each marginal term
      analysed in a model with all other variables), `by = "onedf"` will
      analyse one-degree-of-freedom contrasts sequentially. The argument
      is passed on to
      [`anova.cca`](https://vegandevs.github.io/vegan/reference/anova.cca.html).

  `parallel`

  :   Number of parallel processes or a predefined socket cluster. With
      `parallel = 1` uses ordinary, non-parallel processing. The
      parallel processing is done with parallel package.

  `na.action`

  :   Handling of missing values on the right-hand-side of the formula
      (see [`na.fail`](https://rdrr.io/r/stats/na.fail.html) for
      explanation and alternatives). Missing values are not allowed on
      the left-hand-side. NB, argument `subset` is not implemented.

  `strata`

  :   Groups within which to constrain permutations. The traditional
      non-movable strata are set as Blocks in the
      [permute](https://CRAN.R-project.org/package=permute) package, but
      some more flexible alternatives may be more appropriate.

## Value

An object of class "adonis" (see
[adonis](https://vegandevs.github.io/vegan/reference/vegan-deprecated.html)).

## Details

Samples where one or more predictors are NA are removed.

## Examples

``` r
res <- urt %>%
  perform_adonis(c("plate", "method"), by = "terms")
res
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> adonis2(formula = as.formula(paste("counts_matrix", formula_RHS, sep = " ~ ")), data = metadata, permutations = permutations, by = "terms")
#>           Df SumOfSqs      R2      F Pr(>F)  
#> plate      1    0.663 0.00905 1.9454  0.018 *
#> method     1    0.677 0.00924 1.9863  0.014 *
#> Residual 211   71.913 0.98171                
#> Total    213   73.253 1.00000                
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
