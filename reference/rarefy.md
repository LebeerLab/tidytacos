# Rarefy the samples to a given number of reads

This function performs rarefying. Make sure that all samples contain at
least the minimum number of reads; otherwise, an error might be thrown.

## Usage

``` r
rarefy(ta, n, replace = F)
```

## Arguments

- ta:

  A tidytacos object.

- n:

  Subsample size for rarefying the community.

- replace:

  Whether to replace the read after it has been selected for the
  subsample so it can be sampled again. Default is FALSE.

## Value

A tidytacos object.

## Examples

``` r
# discard samples with less than 1000 reads
urt_1000 <- urt %>%
  add_total_count() %>%
  filter_samples(total_count >= 1000)

# then rarefy to 1000 reads
urt_1000 <- urt_1000 %>% rarefy(1000)
```
