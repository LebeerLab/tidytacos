# Get beta diversity table

`betas()` returns a tidy table with the beta diversity for each
combination of samples.

## Usage

``` r
betas(ta, unique = T, method = "bray", binary = F, ...)
```

## Arguments

- ta:

  A tidytacos object.

- unique:

  A logical scalar. Avoid redundancy by removing all self sample
  comparisons and keep only one of two pairwise comparisons? Default is
  TRUE.

- method:

  The dissimilarity index. See
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  for all options. Default is "bray".

- binary:

  A logical scalar. Perform presence/absence standardization before
  analysis. See
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html).
  Default is FALSE.

- ...:

  Arguments passed on to
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)

  `x`

  :   Community data matrix.

  `diag`

  :   Compute diagonals.

  `upper`

  :   Return only the upper diagonal.

  `na.rm`

  :   Pairwise deletion of missing observations when computing
      dissimilarities (but some dissimilarities may still be `NA`,
      although calculation is handled).

## Details

This function calculates the beta diversity using the
[`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
function of Vegan. It will report one diversity estimate for each
combination of samples.

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
data <-
  create_tidytacos(x, taxa_are_columns = FALSE)

# Report numbers
numbers <- data %>% betas()
```
