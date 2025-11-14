# Add ordination

`add_ord()` adds the first n dimensions of a dimensionality reduction
method performed on a given dissimilarity matrix as new variables to the
sample table of a tidytacos object.

## Usage

``` r
add_ord(ta, distance = "bray", method = "pcoa", dims = 2, binary = FALSE, ...)
```

## Arguments

- ta:

  A tidytacos object.

- distance:

  The distance indices to use, see
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html).

- method:

  The ordination method to use to calculate coordinates. Choices are
  `pcoa`, `tsne`, `umap`.

- dims:

  The amount of dimensions to reduce the dissimilarities to.

- binary:

  Perform presence/absence standardisation before distance computation
  or not.

- ...:

  Additional arguments to pass to the ordination function: either
  [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html),
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) or
  [`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html).

- pseudocount:

  Optional pseudocount to be used in aitchison distance calculation.

## Value

A tidytacos object with the ordination coordinates added.

## Details

This function calculates the dissimilarities between samples followed by
an ordination analysis. It will then add the first n dimensions to the
sample table of a tidytacos object named "ord1", "ord2", ... This
function will also add relative abundances if not present using
[`add_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_rel_abundance.md).

## See also

Other sample-modifiers:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

Other diversity-metrics:
[`add_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_alpha.md),
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md)

## Examples

``` r
# Initiate counts matrix
x <- matrix(
  c(1500, 1300, 280, 356, 456, 678),
  ncol = 3
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2", "sample3")

# Convert to tidytacos object
data <- create_tidytacos(x,
  taxa_are_columns = FALSE
)

# Add pcoa
data <- data %>%
  add_ord()

# The variances of the ordination dimensions can be accessed with
data$ord_variances
#> [1] 1.000000e+00 5.494303e-16 0.000000e+00
```
