# Perform a Sloan neutral model fit on the taxa in the current tidytacos object

This function executes the [ncm](https://rdrr.io/pkg/pctax/man/ncm.html)
function of the pctax package and returns the result.

## Usage

``` r
perform_neutral_model(ta, taxon_name = taxon_id, model = "nls")
```

## Arguments

- ta:

  A tidytacos object.

- taxon_name:

  the name of the column in the taxa table to perform the fit on.

- model:

  the fit method to be used; choice from "nls" or "mle".

## Value

An object of class "ncm_res" (see
[pctax](https://rdrr.io/pkg/pctax/man/ncm.html)).

## Examples

``` r
res <- urt %>% 
  filter_samples(method == "S", location == "NF") %>%
  perform_neutral_model()
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
plot(res)
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
```
