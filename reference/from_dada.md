# DADA2 to a tidytacos object

`from_dada()` returns a tidytacos object given a seqtab and taxa object
from dada2.

## Usage

``` r
from_dada(seqtab, taxa, taxa_are_columns = TRUE)
```

## Arguments

- seqtab:

  Sequence table, output of
  [`dada2::makeSequenceTable()`](https://rdrr.io/pkg/dada2/man/makeSequenceTable.html).

- taxa:

  Taxa table, output of
  [`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html).

- taxa_are_columns:

  A logical scalar. Are the taxa defined in columns?

## Value

A tidytacos object.

## Details

This function will convert two dada2 objects or files into a tidytacos
object.

## See also

Other import-methods:
[`create_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/create_tidytacos.md),
[`from_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/from_phyloseq.md),
[`read_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/read_tidytacos.md)

## Examples

``` r
seqtab <- readRDS(system.file("extdata", "dada2", "seqtab.rds", package = "tidytacos"))
taxa <- readRDS(system.file("extdata", "dada2", "taxa.rds", package = "tidytacos"))
taco <- from_dada(seqtab, taxa)
```
