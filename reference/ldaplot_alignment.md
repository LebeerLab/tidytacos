# Plot LDA topic alignment results

`ldaplot_alignment()` plots the alignment of LDA topics across models.

## Usage

``` r
ldaplot_alignment(lda_alignment, show_product = TRUE, ...)
```

## Arguments

- lda_alignment:

  The output of
  [`align_lda_topics()`](https://lebeerlab.github.io/tidytacos/reference/align_lda_topics.md).

- show_product:

  Whether to show the product or transport alignment.

- ...:

  Arguments passed on to
  [`alto::plot_alignment`](https://rdrr.io/pkg/alto/man/plot_alignment.html)

  `x`

  :   (required) An alignment class object resulting from
      `align_topics`.

  `rect_gap`

  :   (optional) A float describing how much vertical space to put
      between topics within the same model. The units correspond to
      topic masses. Defaults to 0.2.

  `color_by`

  :   (optional) What should the color of topics and weights encode?
      Defaults to 'path'. Other possible arguments are 'coherence',
      'refinement', or 'topic'.

  `model_name_repair_fun`

  :   (optional) How should names be repaired before plotting?

  `label_topics`

  :   (optional, default = `FALSE`) A `logical` specifying if topics
      should be labeled with the `"color_by"` information.

  `add_leaves`

  :   (optional, default = `FALSE`) A `logical` specifying if the topic
      composition of leave-topics should be printed.

  `leaves_text_size`

  :   (optional, default = `10`) specifies the font size of leaves
      annotations in `pt` if `add_leaves` is `TRUE`.

  `n_features_in_leaves`

  :   (optional, default = 3) specifies the maximum number of features
      that should be included in the leaves annotations if `add_leaves`
      is `TRUE`.

  `min_feature_prop`

  :   (optional, default = 0.1) specifies the minimum proportion of a
      feature in a topic for that feature to be included in the leaves
      annotations if `add_leaves` is `TRUE`.

  `top_n_edges`

  :   (optional, `integer`, default = `NULL`) specifies the number of
      edges that should be drawn between the topics of subsequent
      models. The `top_n_edges` with the highest weights are drawn. If
      `NULL` (default), all edges are drawn.

## Value

A ggplot object.
