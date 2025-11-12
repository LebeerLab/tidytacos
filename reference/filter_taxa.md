# Filter the taxa

Filter the taxa

## Usage

``` r
filter_taxa(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Filter criteria for the taxa table.

## Value

A tidytacos object.

## Examples

``` r
# keep only bacterial reads
leaf <- leaf %>% filter_taxa(kingdom == "Bacteria")
```
