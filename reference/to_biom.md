# Write the counts of the tidytacos object to a biom file

Uses taxon_id and sample_id columns to create a dense biom file (v1,
json). By default this is on ASV/OTU level. To do it at any other
taxonomy level, one first needs to aggregate the taxa.

## Usage

``` r
to_biom(ta, filename = "asvs.biom")
```

## Arguments

- ta:

  A tidytacos object.

- filename:

  The name of the resulting biom table file, defaults to 'asvs.biom'.

## See also

Other export-methods:
[`as_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/as_phyloseq.md),
[`to_fasta()`](https://lebeerlab.github.io/tidytacos/reference/to_fasta.md),
[`write_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/write_tidytacos.md)
