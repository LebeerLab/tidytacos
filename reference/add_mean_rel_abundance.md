# Add average relative abundances to taxa table

`add_mean_rel_abundance()` adds mean relative abundance values for each
taxon to the taxa table, overall or per sample group.

## Usage

``` r
add_mean_rel_abundance(ta, condition = NULL, test = NULL)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  A condition variable (character).

- test:

  Differential abundance test to perform.

## Value

A tidytacos object

## Details

If `condition` is specified, the mean relative abundances will be
calculated separately for each group defined by the condition variable.
This variable should be present in the sample table.

If `condition` is specified, differential abundance testing can be
performed by setting the `test` argument. Options are `NULL` (default),
`"wilcox"` or `"t-test"`.

## See also

Other taxa-modifiers:
[`add_prevalence()`](https://lebeerlab.github.io/tidytacos/reference/add_prevalence.md),
[`add_taxon_name()`](https://lebeerlab.github.io/tidytacos/reference/add_taxon_name.md),
[`add_taxon_name_color()`](https://lebeerlab.github.io/tidytacos/reference/add_taxon_name_color.md)
