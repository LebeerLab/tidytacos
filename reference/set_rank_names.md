# Set rank names for a tidytacos object

The rank names are used to describe the taxa in the taxa table.
Tidytacos expects a vector of different names of the ranks used to
describe the taxa. The order of the names should be from the highest
rank to the lowest rank.

## Usage

``` r
set_rank_names(ta, rank_names)
```

## Arguments

- ta:

  a tidytacos object

- rank_names:

  a vector containing the names of the ranks used to describe the taxa

## Value

An updated tidytacos object.

## Details

Eg: c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
