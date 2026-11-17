# Calculate LDA model perplexities for a range of topic numbers

`calculate_lda_perplexities()` estimates LDA models for a range of topic
numbers and returns the perplexity of each model.

## Usage

``` r
calculate_lda_perplexities(ta, k_range = 2:10, control = list(seed = 42), ...)
```

## Arguments

- ta:

  A tidytacos object.

- k_range:

  A range of topic numbers to estimate.

- control:

  An optional list of arguments to pass to the
  [`topicmodels::LDA()`](https://rdrr.io/pkg/topicmodels/man/lda.html)
  function, default includes a seed set at 42

- ...:

  Arguments passed on to
  [`perform_lda`](https://lebeerlab.github.io/tidytacos/reference/perform_lda.md)

  `k`

  :   The number of topics to estimate.

  `min_prevalence`

  :   The lowest percentage (0-1) of samples taxa need to be present in
      for the taxa to be used in model estimation.

  `taxon`

  :   The column name in the taxa table with taxa identifiers.

  `sample`

  :   The column name in the sample table with sample identifiers.

- seed:

  A seed for reproducibility.

## Value

A tibble with the perplexity of each model and the number of topics.
