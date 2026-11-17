# LDA model estimation

`perform_lda()` estimates a Latent Dirichlet Allocation (LDA) model on
the counts matrix of a tidytacos object. The function returns the
estimated topics and terms, as well as the perplexity of the model.

## Usage

``` r
perform_lda(
  ta,
  k,
  min_prevalence = 0.05,
  taxon = taxon_id,
  sample = sample_id,
  control = list(seed = 42),
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- k:

  The number of topics to estimate.

- min_prevalence:

  The lowest percentage (0-1) of samples taxa need to be present in for
  the taxa to be used in model estimation.

- taxon:

  The column name in the taxa table with taxa identifiers.

- sample:

  The column name in the sample table with sample identifiers.

- control:

  An optional list of arguments to pass to the
  [`topicmodels::LDA()`](https://rdrr.io/pkg/topicmodels/man/lda.html)
  function, default includes a seed set at 42

- ...:

  Arguments passed on to
  [`topicmodels::LDA`](https://rdrr.io/pkg/topicmodels/man/lda.html)

  `x`

  :   Object of class `"DocumentTermMatrix"` with term-frequency
      weighting or an object coercible to a `"simple_triplet_matrix"`
      with integer entries.

  `method`

  :   The method to be used for fitting; currently `method = "VEM"` or
      `method= "Gibbs"` are supported.

  `model`

  :   Object of class `"LDA"` for initialization.

## Value

A list of estimated topics, terms and the perplexity of the model.
