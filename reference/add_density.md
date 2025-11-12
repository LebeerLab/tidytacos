# Add density to count table

`add_density()` adds densities (bacterial biomass per sample mass or
volume) to the count table of a tidytacos object under the column name
"density". Can only be used of a taxon was spiked into samples during
library prep.

## Usage

``` r
add_density(
  ta,
  spike_taxon,
  spike_added = spike_added,
  material_sampled = material_sampled
)
```

## Arguments

- ta:

  A tidytacos object.

- spike_taxon:

  The taxon id of the spike.

- spike_added:

  The column name of the sample table which indicates how much spike was
  added per sample, e.g. 16S rRNA gene copy numbers added to the DNA
  extraction tube.

- material_sampled:

  The column name indicating the amount of material from which DNA was
  extracted, e.g gram of soil. This parameter encourages researchers to
  consider that absolute abundances are only meaningful if they can be
  translated into densities.

## Value

A tidytacos object with densities added to the count table.

## See also

Other count-modifiers:
[`add_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_absolute_abundance.md),
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
data$samples$grams_source <- c(4, 5)

# Add density
data <- data %>%
  add_density(spike_taxon = "t3", material_sampled=grams_source)
```
