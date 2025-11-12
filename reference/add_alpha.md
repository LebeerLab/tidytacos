# Add alpha diversity measure

`add_alpha()` adds an alpha diversity measures to the sample table of a
tidytacos object.

## Usage

``` r
add_alpha(
  ta,
  method = "invsimpson",
  keep_empty_samples = FALSE,
  subsample = FALSE,
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- method:

  The diversity measure to use, see
  [`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html)
  for further information on these.

- keep_empty_samples:

  Whether to keep empty samples or not.

- subsample:

  Wether to use subsampling to be able to compare samples of varying
  sequencing depths. [Schloss et al.,
  2023](https://journals.asm.org/doi/10.1128/msphere.00355-23)

- ...:

  Arguments passed on to
  [`add_subsampled_alpha`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md)

  `min_lib_size`

  :   the minimum lib size samples need to have. Samples with lower lib
      sizes will be discarded and samples with a higher readcount will
      be itteratively subsampled to this readcount to allow for a fair
      comparison across read_depths.

  `itterations`

  :   the amount of itterations for subsampling. Please report this
      number in your research.

## Value

A tidytacos object with the alpha diversity measure added.

## Details

This function can add different alpha diversity measures to the sample
table, specified by the method argument. The following methods are
available:

- invsimpson: Inverse Simpson index

- shannon: Shannon index

- simpson: Simpson index

- pielou: Pielou's evenness index

- obs: Observed richness

- s.chao1: Chao1 richness estimator

- s.ace: ACE richness estimator

## See also

Other sample-modifiers:
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_dominant_taxa()`](https://lebeerlab.github.io/tidytacos/reference/add_dominant_taxa.md),
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_sample_clustered()`](https://lebeerlab.github.io/tidytacos/reference/add_sample_clustered.md),
[`add_spike_ratio()`](https://lebeerlab.github.io/tidytacos/reference/add_spike_ratio.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
[`add_total_absolute_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_total_absolute_abundance.md),
[`add_total_count()`](https://lebeerlab.github.io/tidytacos/reference/add_total_count.md),
[`add_total_density()`](https://lebeerlab.github.io/tidytacos/reference/add_total_density.md),
[`cluster_samples()`](https://lebeerlab.github.io/tidytacos/reference/cluster_samples.md)

Other diversity-metrics:
[`add_alphas()`](https://lebeerlab.github.io/tidytacos/reference/add_alphas.md),
[`add_ord()`](https://lebeerlab.github.io/tidytacos/reference/add_ord.md),
[`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md)

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

# Add alpha diversity measures
data <- data %>%
  add_alpha()

data <- data %>%
  add_alpha(method = "shannon")
```
