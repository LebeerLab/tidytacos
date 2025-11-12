# Plot LDA beta values

`ldaplot_beta()` plots the beta values of LDA topics.

## Usage

``` r
ldaplot_beta(lda_alignment, show_product = TRUE, ...)
```

## Arguments

- lda_alignment:

  The output of
  [`align_lda_topics()`](https://lebeerlab.github.io/tidytacos/reference/align_lda_topics.md).

- show_product:

  Whether to show the product or transport alignment.

- ...:

  Arguments passed on to
  [`alto::plot_beta`](https://rdrr.io/pkg/alto/man/plot_beta.html)

  `x`

  :   (required) An alignment class object resulting from
      `align_topics`.

  `models`

  :   Which models to display in the heatmap? Defaults to `"all"`,
      meaning that all models are shown. If given `"last"`, only the
      last model in the models list will be plotted. If given a vector
      of characters, it will plot only models whose names in the
      original models list match. Similarly, if given a list of
      integers, only the models lying at those indices in the original
      model list will be visualized.

  `filter_by`

  :   (optional, default = `"beta"`) a character specifying if the data
      (beta matrices) should be filtered by the average `"beta"` across
      topics or by the `"distinctiveness"` of the features.

  `x_axis`

  :   (optional, default = `"index"`) a character specifying if the
      x-axis should display topic indices (`"index"`) such that they
      match the alignment plot order or topic names (`"label"`).

  `threshold`

  :   (optional, default = 0.001) Words (features) with less than this
      average beta or distinctiveness across all topics are ignored

  `n_features`

  :   (optional) alternative to `threshold`. The maximum number of words
      (features) to display along rows of the plot.

  `beta_aes`

  :   Should word probabilities within a topic be encoded using circle
      size (`"size"`) or opacity (`"alpha"`) ? Defaults to `"size"`.

  `color_by`

  :   (optional) What should the color of topics and weights encode?
      Defaults to 'path'. Other possible arguments are 'coherence',
      'refinement', or 'topic'.

## Value

A ggplot object.
