# Aligns LDA topics across models

`align_lda_topics()` estimates LDA models for a range of topic numbers
and aligns the topics across models using the product and transport
method.

## Usage

``` r
align_lda_topics(
  ta,
  k_range = 2:10,
  seed = 42,
  min_prevalence = 0.05,
  taxon = taxon_id,
  sample = sample_id,
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- k_range:

  A range of topic numbers to estimate.

- seed:

  A seed for reproducibility.

- min_prevalence:

  The lowest percentage (0-1) of samples taxa need to be present in for
  the taxa to be used in model estimation.

- taxon:

  The column name in the taxa table with taxa identifiers.

- sample:

  The column name in the sample table with sample identifiers.

- ...:

  Arguments passed on to
  [`alto::run_lda_models`](https://rdrr.io/pkg/alto/man/run_lda_models.html)

  `data`

  :   (required) a `matrix`, `data.frame` or
      [`slam::simple_triplet_matrix`](https://rdrr.io/pkg/slam/man/matrix.html)
      containing the counts (integers) of each feature (e.g. words) and
      each sample (or document). If data is provided as `matrix` or
      `data.frame`, each row is a sample, each column is a feature.

  `lda_varying_params_lists`

  :   (required) a `list` specifying the parameter for each models that
      needs to be ran. Currently, supported parameters are "k" (the
      number of topic), "method" ("VEM" or "Gibbs"), and "control", a
      list of type `LDAcontrol`. See
      [`topicmodels::LDA`](https://rdrr.io/pkg/topicmodels/man/lda.html)
      for details and below for examples.

  `lda_fixed_params_list`

  :   (optional) a `list` specifying the parameters common to all models
      to be fitted. Values provided by `lda_fixed_params_list` are
      overwritten by those provided by `lda_varying_params_lists`.

  `dir`

  :   (optional) a `character` specifying the directory in which
      individual LDA models should be stored. If not specified,
      individual LDA models are not stored. This option is especially
      useful for data exploration as it allows to save execution time if
      one wishes to add models to an existing model list. (see examples)

  `reset`

  :   (optional, default = `FALSE`). Should any cached models in the
      save directory be cleared?

  `verbose`

  :   (optional, default = `FALSE`) Print verbose output while running
      models?

## Value

A list of estimated models, aligned topics using the product method and
aligned topics using the transport method.
