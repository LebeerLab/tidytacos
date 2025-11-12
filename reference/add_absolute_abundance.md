# Add absolute abundances to count table

`add_absolute_abundance()` calculates absolute abundances of taxa in
samples given a taxon that was spiked into all of the samples during
library prep. The function then adds these absolute abundances to the
count table of the tidytacos object under the column name
"absolute_abundance".

## Usage

``` r
add_absolute_abundance(ta, spike_taxon, spike_added = spike_added)
```

## Arguments

- ta:

  A tidytacos object.

- spike_taxon:

  The taxon id of the spike.

- spike_added:

  The column name in the sample table which indicates how much spike was
  added per sample.

## Value

A tidytacos object with absolute abundances added to the count table.

## See also

Other count-modifiers:
[`add_density()`](https://lebeerlab.github.io/tidytacos/reference/add_density.md),
[`add_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_rel_abundance.md)

## Examples

``` r
# Initiate count matrix
x <- matrix(
  c(1500, 1300, 14, 280, 356, 9),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2", "taxon3")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x,
  taxa_are_columns = FALSE
)
data$samples$spike_added <- c(100, 150)

# Add total abundance
data <- data %>%
  add_absolute_abundance(spike_taxon = "t3")
```
