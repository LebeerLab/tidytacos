# Filter the samples

Filter the samples

## Usage

``` r
filter_samples(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Filter criteria for the samples table.

## Value

A tidytacos object.

## Examples

``` r
# subset urt to keep only nasopharynx samples
urt_nf <- urt %>% filter_samples(location == "NF")
# subset urt to keep only samples from plate 1 and 2
urt_plate_1_2 <- urt %>% filter_samples(plate %in% c(1, 2))
# subset the blanks in leaf
leaf_blanks <- leaf %>% filter_samples(startsWith(description, "BLANK"))
```
