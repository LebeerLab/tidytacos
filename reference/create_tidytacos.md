# Initiate tidytacos object

`create_tidytacos()` returns a tidytacos object given a numeric matrix.

## Usage

``` r
create_tidytacos(
  counts_matrix,
  taxa_are_columns = TRUE,
  allow_non_count = FALSE
)
```

## Arguments

- counts_matrix:

  Numerical matrix containing the count data.

- taxa_are_columns:

  A logical scalar. Are the taxa defined in columns?

- allow_non_count:

  Allow counts that are less than or equal to 0.

## Value

A tidytacos object.

## Details

This function initiates a tidytacos object based on a numeric matrix. It
will automatically create a dummy taxa table and sample table which will
need to be updated using the function
[`add_metadata()`](https://lebeerlab.github.io/tidytacos/reference/add_metadata.md).
When the taxa table is updated, the rank names can be set using the
function
[`set_rank_names()`](https://lebeerlab.github.io/tidytacos/reference/set_rank_names.md)
to make the tidytacos object aware of the taxonomy level order (from
high to low). The taxa table should contain at the very least one rank
name. The default rank names used by tidytacos are "domain", "phylum",
"class", "order", "family", "genus" and "species".

## See also

Other import-methods:
[`from_dada()`](https://lebeerlab.github.io/tidytacos/reference/from_dada.md),
[`from_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/from_phyloseq.md),
[`read_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/read_tidytacos.md)

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
# Add taxonomy information
taxonomy <- tibble::tibble(
  taxon = c("taxon1", "taxon2"),
  domain = c("Bacteria", "Bacteria"),
  phylum = c("Bacillota", "Pseudomonadota"),
  class = c("Bacilli", "Gammaproteobacteria"),
  order = c("Lactobacillales", "Enterobacteriales"),
  family = c("Lactobacillaceae", "Enterobacteriaceae"),
  genus = c("Lactobacillus", "Escherichia"),
  species = c("Lactobacillus crispatus", "Escherichia coli")
)

data <- add_metadata(data, taxonomy, table_type = "taxa")
#> Joining with `by = join_by(taxon)`
# rank names are inferred from the taxa table, but can be set manually
# if we're not happy with the inferred rank names.
data <- set_rank_names(
  data, c("domain", "phylum", "class", "order", "family", "genus")
)
```
