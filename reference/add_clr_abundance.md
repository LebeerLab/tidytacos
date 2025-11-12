# Perform a centered log ratio transformation on the readcounts.

`add_clr_abundance()` calculates the log ration transformed values for
each taxon in each sample and adds these data in a new table,
clr_counts. Alternatively, using 'overwrite', the clr transformed data
can replace the 'counts' column in the count table.

## Usage

``` r
add_clr_abundance(ta, overwrite = F, pseudocount = 1, ...)
```

## Arguments

- ta:

  A tidytacos object.

- overwrite:

  Whether or not the counts table is to be overwritten with the
  transformed counts.

- pseudocount:

  A pseudocount to be added to the counts before transformation. If
  false or zero will perform robust CLR.

- ...:

  Arguments passed on to
  [`vegan::decostand`](https://vegandevs.github.io/vegan/reference/decostand.html)

  `x`

  :   Community data, a matrix-like object. For `decobackstand`
      standardized data.

  `method`

  :   Standardization method. See Details for available options.

  `MARGIN`

  :   Margin, if default is not acceptable. `1` = rows, and `2` =
      columns of `x`.

  `range.global`

  :   Matrix from which the range is found in `method = "range"`. This
      allows using same ranges across subsets of data. The dimensions of
      `MARGIN` must match with `x`.

  `logbase`

  :   The logarithm base used in `method = "log"`.

  `na.rm`

  :   Ignore missing values in row or column standardizations. The `NA`
      values remain as `NA`, but they are ignored in standardization of
      other values.

## Value

A tidytacos object.
