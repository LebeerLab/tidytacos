# Create extra variables in the sample table

Create extra variables in the sample table

## Usage

``` r
mutate_samples(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Mutate criteria for the samples table.

## Value

A tidytacos object.

## Examples

``` r
# change the sample column to lowercase
urt <- urt %>% mutate_samples(sample = tolower(sample))
```
