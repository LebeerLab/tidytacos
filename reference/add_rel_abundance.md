# Add relative abundances to count table

`add_rel_abundance()` calculates relative abundances of taxa in samples
and adds them to the count table of a tidytacos object under the column
name "rel_abundance".

## Usage

``` r
add_rel_abundance(ta)
```

## Arguments

- ta:

  A tidytacos object.

## Value

A tidytacos object with relative abundances added to the count table.

## See also

Other count-modifiers:
[`add_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_absolute_abundance.md),
[`add_density()`](https://lebeerlab.github.io/tidytacos/reference/add_density.md)

## Examples

``` r
# Initiate count matrix
x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x,
  taxa_are_columns = FALSE
)

# Add relative abundance
data <- data %>% add_rel_abundance()
```
