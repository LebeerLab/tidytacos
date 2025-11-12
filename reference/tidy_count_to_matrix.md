# Convert counts tidy data frame to matrix

`tidy_count_to_matrix()` returns a numerical matrix given a tidy counts
data frame.

## Usage

``` r
tidy_count_to_matrix(counts, value = count)
```

## Arguments

- counts:

  The counts tidy data frame that will be converted.

- value:

  Name of column containing the counts data. Default is "counts".

## Details

This function will convert a counts tidy data frame into a numerical
counts matrix. To convert a numerical counts matrix into a counts tidy
data frame use `tidy_count_to_matrix()`.
