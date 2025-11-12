# Add logratios

`add_logratio()` computes pairwise logratio values between all taxa and
adds these to the tidytacos object in the form of a table called
logratios.

## Usage

``` r
add_logratio(ta, max_taxa = 50)
```

## Arguments

- ta:

  A tidytacos object.

- max_taxa:

  The maximum number of taxa to use.

## Value

A tidytacos object with an extra table logratios

## Details

If `max_taxa` is smaller than the number of taxa in the dataset, the
taxa with the highest prevalence will be selected.

IMPORTANT: this function adds pseudocounts of one to all abundances
before calculating the logratios.
