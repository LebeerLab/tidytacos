# Write community data in tidytacos format

`write_tidytacos()` saves the tidytacos object into 3 .csv files. This
format allows easy loading of the tidytacos object using the
[`read_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/read_tidytacos.md)
function.

## Usage

``` r
write_tidytacos(ta, dout)
```

## Arguments

- ta:

  A tidytacos object.

- dout:

  The directory to store the three tidytacos tables in.

## See also

Other export-methods:
[`as_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/as_phyloseq.md),
[`to_biom()`](https://lebeerlab.github.io/tidytacos/reference/to_biom.md),
[`to_fasta()`](https://lebeerlab.github.io/tidytacos/reference/to_fasta.md)
