# Convert tidytacos object to phyloseq object

`as_phyloseq()` returns a phyloseq object given a tidytacos object.

## Usage

``` r
as_phyloseq(ta, sample = sample, taxon = taxon_id)
```

## Arguments

- ta:

  A tidytacos object.

- sample:

  The sample names required for a phyloseq object. Default is "sample"
  column of the sample table of the tidytacos object.

- taxon:

  The taxon names required for a phyloseq object. Default is the
  "taxon_id" column in the taxon table of the tidytacos object.

## Details

This function will convert a tidytacos object into a phyloseq object for
alternative processing using the phyloseq package. To convert from a
phyloseq object to a tidytacos object use
[`from_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/from_phyloseq.md).

## See also

Other export-methods:
[`to_biom()`](https://lebeerlab.github.io/tidytacos/reference/to_biom.md),
[`to_fasta()`](https://lebeerlab.github.io/tidytacos/reference/to_fasta.md),
[`write_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/write_tidytacos.md)
