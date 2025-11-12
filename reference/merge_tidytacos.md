# Merge two tidytacos objects

`merge_tidytacos()` merges two tidytacos objects and returns one single
tidytacos object.

## Usage

``` r
merge_tidytacos(ta1, ta2, taxon_identifier = sequence)
```

## Arguments

- ta1:

  The first tidytacos object.

- ta2:

  The second tidytacos object.

- taxon_identifier:

  The column name in the taxa tables which identify unique taxa. Default
  is sequence.

## Details

This function will merge two tidytacos objects into one. It is useful if
one wants to merge data obtained from different sequencing runs.
Therefore, this function requires that both tidytacos objects contain a
"run" variable in their samples table, indicating their origin.
