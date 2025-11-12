# Filter the counts

Filter the counts

## Usage

``` r
filter_counts(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Filter criteria for the counts table.

## Value

A tidytacos object.

## Examples

``` r
# remove singletons
urt <- urt %>% filter_counts(count > 1)
```
