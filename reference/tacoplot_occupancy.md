# Return an occupancy-abundance plot to visualize the core taxa of a community.

Return an occupancy-abundance plot to visualize the core taxa of a
community.

## Usage

``` r
tacoplot_occupancy(ta, taxon_name = taxon_id, ...)
```

## Arguments

- ta:

  A tidytacos object of the community in question.

- taxon_name:

  the name of the taxa table column to perform the fit on.

- ...:

  Arguments passed on to
  [`perform_neutral_model`](https://lebeerlab.github.io/tidytacos/reference/perform_neutral_model.md)

  `model`

  :   the fit method to be used; choice from "nls" or "mle".

## Examples

``` r
urt %>% tacoplot_occupancy()
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Joining with `by = join_by(taxon_id)`
```
