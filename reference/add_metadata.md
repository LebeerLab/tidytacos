# Add metadata to the tidytacos object

`add_metadata()` adds sample or taxon metadata to the sample or taxon
table, respectively, of a tidytacos object.

## Usage

``` r
add_metadata(ta, metadata, table_type = "sample")
```

## Arguments

- ta:

  A tidytacos object.

- metadata:

  A tibble containing data for each sample or taxon. Samples/taxa should
  be rows, while metadata variables should be columns. At least one
  column name needs to be shared with the sample or taxa table of the
  tidytacos object. The default shared column name is 'sample' for
  samples and 'taxon' for taxa.

- table_type:

  The type of table to add, either 'sample' or 'taxa'.

## Value

A tidytacos object with the metadata added.

A tidytacos object with metadata columns added to the taxa or sample
table.

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

## Examples

``` r
# Initiate counts matrix
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

# Initiate sample tibble
sample <- c("sample1", "sample2")
environment <- c("food fermentation", "human stool")
sample_tibble <- tibble::tibble(sample, environment)

# Add sample tibble to tidytacos object
data <- data %>%
  add_metadata(sample_tibble)
#> Joining with `by = join_by(sample)`

# Initiate taxon tibble
genera <- c("Lactobacillus", "Limosilactobacillus")
species <- c("crispatus", "reuteri")
taxonomy <- tibble::tibble(taxon = rownames(x), genera, species)

# Add taxon tibble to tidytacos object
data <- data %>%
  add_metadata(taxonomy, table_type = "taxa")
#> Joining with `by = join_by(taxon)`
```
