# Read community data written by tidytacos

`read_tidytacos()` reads the three .csv files created by the
[`write_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/write_tidytacos.md)
function and returns a tidytacos object.

## Usage

``` r
read_tidytacos(
  din,
  samples = "samples.csv",
  taxa = "taxa.csv",
  counts = "counts.csv"
)
```

## Arguments

- din:

  directory containing the sample, taxa and counts table in csv format

- samples:

  the name of the samples table, defaults to samples.csv

- taxa:

  the name of the taxa table, defaults to taxa.csv

- counts:

  the name of the counts table, defaults to counts.csv

## Value

A tidytacos object.

## Details

The samples.csv file should contain a column named "sample_id". The
taxa.csv file should contain a column named "taxon_id" and at the very
least one rank name. The default rank names used by tidytacos are
"domain", "phylum", "class", "order", "family", "genus" and "species".
The counts.csv file should contain columns named "sample_id", "taxon_id"
and "count".

## See also

Other import-methods:
[`create_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/create_tidytacos.md),
[`from_dada()`](https://lebeerlab.github.io/tidytacos/reference/from_dada.md),
[`from_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/from_phyloseq.md)
