# Perform network inference with SparCC on a tidytacos object, after dropping rare taxa. See [`SpiecEasi::sparcc()`](https://rdrr.io/pkg/SpiecEasi/man/sparcc.html).

Perform network inference with SparCC on a tidytacos object, after
dropping rare taxa. See
[`SpiecEasi::sparcc()`](https://rdrr.io/pkg/SpiecEasi/man/sparcc.html).

## Usage

``` r
network(
  ta,
  min_occurrence = 0.05,
  taxon_name = taxon,
  sample_name = sample,
  calculate_p = FALSE,
  ...
)
```

## Arguments

- ta:

  a tidytacos object

- min_occurrence:

  Percentage of samples the taxon needs to be present in for it to be
  considered in the analysis.

- taxon_name:

  Column name of the taxon identifier, by default taxon.

- sample_name:

  Column name of the sample identifier, by default sample. considered
  zero by the inner SparCC loop.

- calculate_p:

  whether to calculate p-values or not. This can be time consuming due
  to the many iterations needed. Iterations can be set with the R
  parameter and multiple cores through ncpus.

- ...:

  Arguments passed on to
  [`SpiecEasi::sparcc`](https://rdrr.io/pkg/SpiecEasi/man/sparcc.html)

  `data`

  :   Community count data matrix

  `iter`

  :   Number of iterations in the outer loop

  `inner_iter`

  :   Number of iterations in the inner loop

  `th`

  :   absolute value of correlations below this threshold are considered
      zero by the inner SparCC loop.
