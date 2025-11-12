# Determine the correlation between the distance of the counts in a tidytacos object and a sample variable, multiple sample variables or another matrix.

This function performs a mantel test using the dissimilarity matrix of
the tidytacos object supplied and a second distance matrix generated
from user input.

## Usage

``` r
perform_mantel_test(ta, comparison, ...)
```

## Arguments

- ta:

  A tidytacos object.

- comparison:

  A distance to compare against. This can be any of the following:

  - The name of the variable in the sample table to use for comparison

  - A list of names of variables in the sample table.

  - A distance matrix.

- ...:

  Additional arguments to pass to the mantel function.

## Value

The mantel test statistics
