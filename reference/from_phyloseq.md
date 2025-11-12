# Convert phyloseq object to tidytacos object

`from_phyloseq()` returns a tidytacos object given a phyloseq object.

## Usage

``` r
from_phyloseq(ps)
```

## Arguments

- ps:

  Phyloseq object.

## Value

A tidytacos object.

## Details

This function will convert a phyloseq object into a tidytacos object. To
convert from a tidytacos object to a phyloseq object use
[`as_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/as_phyloseq.md).

## See also

Other import-methods:
[`create_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/create_tidytacos.md),
[`from_dada()`](https://lebeerlab.github.io/tidytacos/reference/from_dada.md),
[`read_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/read_tidytacos.md)

## Examples

``` r
phylo_obj <- readRDS(system.file("extdata","phyloseq.rds",package='tidytacos'))
taco <- from_phyloseq(phylo_obj)
#> Joining with `by = join_by(sample)`
#> Joining with `by = join_by(taxon)`
```
