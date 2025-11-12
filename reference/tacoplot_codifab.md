# Generate a compositional differential abundance plot

This function returns a plot to visualize differential abundance of taxa
between conditions, compared to all other taxa as references. These
differential abundances should already have been calculated with
[`add_codifab()`](https://lebeerlab.github.io/tidytacos/reference/add_codifab.md).
Taxa that have a relatively high number of significantly different
ratios, can be considered more abundant in one condition versus the
other.

## Usage

``` r
tacoplot_codifab(ta, diffabun_var)
```

## Arguments

- ta:

  A tidytacos object.

- diffabun_var:

  The variable with differential abundances in the taxon_pair table.

## Value

A ggplot object

## Details

Significance of tests is determined by capping the false discovery rate
at 10%, using the method of Benjamini and Yekutieli, which is developed
for non-independent tests. See
[p.adjust](https://rdrr.io/r/stats/p.adjust.html).

## See also

Other codifab-functions:
[`add_codifab()`](https://lebeerlab.github.io/tidytacos/reference/add_codifab.md)
