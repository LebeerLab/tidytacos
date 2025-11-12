# Return a boxplot of every alpha metric per group in the samples table of a tidytaco object. If no alpha metrics are present, all available ones are added.

Return a boxplot of every alpha metric per group in the samples table of
a tidytaco object. If no alpha metrics are present, all available ones
are added.

## Usage

``` r
tacoplot_alphas(
  ta,
  group_by,
  compare_means = FALSE,
  keep_empty_samples = FALSE,
  subsample_metric = NULL,
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- group_by:

  The name of a variable in the samples table on which to group the
  samples.

- compare_means:

  Add the result of a statistical test to the plot, comparing the means
  of the groups. Default is FALSE.

- keep_empty_samples:

  Whether to discard samples not containing any counts or not. By
  default these are removed.

- subsample_metric:

  if the alpha diversities are precalculated with
  [`add_subsampled_alpha()`](https://lebeerlab.github.io/tidytacos/reference/add_subsampled_alpha.md),
  choose here "mean" or "median" to represent the alpha diversity.

- ...:

  Arguments passed on to
  [`ggpubr::stat_compare_means`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html)

  `mapping`

  :   Set of aesthetic mappings created by
      [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html). If
      specified and `inherit.aes = TRUE` (the default), it is combined
      with the default mapping at the top level of the plot. You must
      supply `mapping` if there is no plot mapping.

  `data`

  :   The data to be displayed in this layer. There are three options:

      If `NULL`, the default, the data is inherited from the plot data
      as specified in the call to
      [`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

      A `data.frame`, or other object, will override the plot data. All
      objects will be fortified to produce a data frame. See
      [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
      for which variables will be created.

      A `function` will be called with a single argument, the plot data.
      The return value must be a `data.frame`, and will be used as the
      layer data. A `function` can be created from a `formula` (e.g.
      `~ head(.x, 10)`).

  `method`

  :   a character string indicating which method to be used for
      comparing means.

  `paired`

  :   a logical indicating whether you want a paired test. Used only in
      [`t.test`](https://rdrr.io/r/stats/t.test.html) and in
      [wilcox.test](https://rdrr.io/r/stats/wilcox.test.html).

  `method.args`

  :   a list of additional arguments used for the test method. For
      example one might use
      `method.args = list(alternative = "greater")` for wilcoxon test.

  `ref.group`

  :   a character string specifying the reference group. If specified,
      for a given grouping variable, each of the group levels will be
      compared to the reference group (i.e. control group).

      `ref.group` can be also `".all."`. In this case, each of the
      grouping variable levels is compared to all (i.e. basemean).

  `comparisons`

  :   A list of length-2 vectors. The entries in the vector are either
      the names of 2 values on the x-axis or the 2 integers that
      correspond to the index of the groups of interest, to be compared.

  `hide.ns`

  :   logical value. If TRUE, hide ns symbol when displaying
      significance levels.

  `label.sep`

  :   a character string to separate the terms. Default is ", ", to
      separate the correlation coefficient and the p.value.

  `label`

  :   character string specifying label type. Allowed values include
      "p.signif" (shows the significance levels), "p.format" (shows the
      formatted p value).

  `label.x.npc,label.y.npc`

  :   can be `numeric` or `character` vector of the same length as the
      number of groups and/or panels. If too short they will be
      recycled.

      - If `numeric`, value should be between 0 and 1. Coordinates to be
        used for positioning the label, expressed in "normalized parent
        coordinates".

      - If `character`, allowed values include: i) one of c('right',
        'left', 'center', 'centre', 'middle') for x-axis; ii) and one of
        c( 'bottom', 'top', 'center', 'centre', 'middle') for y-axis.

  `label.x,label.y`

  :   `numeric` Coordinates (in data units) to be used for absolute
      positioning of the label. If too short they will be recycled.

  `vjust`

  :   move the text up or down relative to the bracket.

  `tip.length`

  :   numeric vector with the fraction of total height that the bar goes
      down to indicate the precise column. Default is 0.03. Can be of
      same length as the number of comparisons to adjust specifically
      the tip lenth of each comparison. For example tip.length = c(0.01,
      0.03).

      If too short they will be recycled.

  `bracket.size`

  :   Width of the lines of the bracket.

  `step.increase`

  :   numeric vector with the increase in fraction of total height for
      every additional comparison to minimize overlap.

  `symnum.args`

  :   a list of arguments to pass to the function
      [`symnum`](https://rdrr.io/r/stats/symnum.html) for symbolic
      number coding of p-values. For example,
      `symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))`.

      In other words, we use the following convention for symbols
      indicating statistical significance:

      - `ns`: p \> 0.05

      - `*`: p \<= 0.05

      - `**`: p \<= 0.01

      - `***`: p \<= 0.001

      - `****`: p \<= 0.0001

  `geom`

  :   The geometric object to use to display the data, either as a
      `ggproto` `Geom` subclass or as a string naming the geom stripped
      of the `geom_` prefix (e.g. `"point"` rather than `"geom_point"`)

  `position`

  :   Position adjustment, either as a string naming the adjustment
      (e.g. `"jitter"` to use `position_jitter`), or the result of a
      call to a position adjustment function. Use the latter if you need
      to change the settings of the adjustment.

  `na.rm`

  :   If FALSE (the default), removes missing values with a warning. If
      TRUE silently removes missing values.

  `show.legend`

  :   logical. Should this layer be included in the legends? `NA`, the
      default, includes if any aesthetics are mapped. `FALSE` never
      includes, and `TRUE` always includes. It can also be a named
      logical vector to finely select the aesthetics to display.

  `inherit.aes`

  :   If `FALSE`, overrides the default aesthetics, rather than
      combining with them. This is most useful for helper functions that
      define both data and aesthetics and shouldn't inherit behaviour
      from the default plot specification, e.g.
      [`borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).
