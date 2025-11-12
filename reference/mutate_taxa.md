# Create extra variables in the taxa table

Create extra variables in the taxa table

## Usage

``` r
mutate_taxa(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Mutate criteria for the taxa table.

## Value

A tidytacos object.

## Examples

``` r
urt <- urt %>% mutate_taxa(species = paste(genus, species))
```
