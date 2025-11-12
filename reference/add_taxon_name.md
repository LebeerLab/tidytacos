# Add sensible taxon name to taxon table

`add_taxon_name()` creates sensible taxa names by -under default
conditions- combining the genus name with a number. The number is only
added if there is more than one taxon of that genus. The number
indicates the rank of abundance, with 1 indicating the taxon with the
highest mean relative abundance within the genus. If genus
classification is not available the next most detailed taxonomic rank
which is available is used. The sensible taxon name is added to the
taxon table under the column name `taxon_name`.

## Usage

``` r
add_taxon_name(ta, method = "mean_rel_abundance", include_species = FALSE)
```

## Arguments

- ta:

  A tidytacos object.

- method:

  The method by which to arrange the taxon names. Currently only
  mean_rel_abundance.

- include_species:

  Whether to include the species name or not.

## Value

A tidytacos object.

## See also

Other taxa-modifiers:
[`add_mean_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_mean_rel_abundance.md),
[`add_prevalence()`](https://lebeerlab.github.io/tidytacos/reference/add_prevalence.md),
[`add_taxon_name_color()`](https://lebeerlab.github.io/tidytacos/reference/add_taxon_name_color.md)

## Examples

``` r
urt_g <- urt %>% add_taxon_name()
# add the species name if present (which is often uncertain in amplicon data)
urt_s <- urt %>% add_taxon_name(include_species = TRUE)
```
