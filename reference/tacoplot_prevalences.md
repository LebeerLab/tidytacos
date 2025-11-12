# Return a heatmap of prevalence of taxa in groups of samples

Return a heatmap of all taxa above a certain threshold prevalence per
condition, clusters them and compares prevalences with a fisher_test.

## Usage

``` r
tacoplot_prevalences(
  ta,
  condition,
  cutoff = 0.1,
  fisher = T,
  adjp_method = "fdr",
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  The row name of the condition which rrevalences are to be compared.

- cutoff:

  The minimum prevalence of a taxon to be included in the heatmap.

- fisher:

  Run a fisher test on the relative prevalences in each condition and
  plot the resulting adjusted p-values as \*(\<.05), \*\*(\<.01),
  \*\*\*(\<.001) or \*\*\*\*(\<.0001).

- adjp_method:

  The method to adjust the p-values, see
  [`rstatix::adjust_pvalue()`](https://rpkgs.datanovia.com/rstatix/reference/adjust_pvalue.html).

- ...:

  Arguments passed on to
  [`pheatmap::pheatmap`](https://rdrr.io/pkg/pheatmap/man/pheatmap.html)

  `mat`

  :   numeric matrix of the values to be plotted.

  `color`

  :   vector of colors used in heatmap.

  `kmeans_k`

  :   the number of kmeans clusters to make, if we want to aggregate the
      rows before drawing heatmap. If NA then the rows are not
      aggregated.

  `breaks`

  :   a sequence of numbers that covers the range of values in mat and
      is one element longer than color vector. Used for mapping values
      to colors. Useful, if needed to map certain values to certain
      colors, to certain values. If value is NA then the breaks are
      calculated automatically. When breaks do not cover the range of
      values, then any value larger than `max(breaks)` will have the
      largest color and any value lower than` min(breaks)` will get the
      lowest color.

  `border_color`

  :   color of cell borders on heatmap, use NA if no border should be
      drawn.

  `cellwidth`

  :   individual cell width in points. If left as NA, then the values
      depend on the size of plotting window.

  `cellheight`

  :   individual cell height in points. If left as NA, then the values
      depend on the size of plotting window.

  `scale`

  :   character indicating if the values should be centered and scaled
      in either the row direction or the column direction, or none.
      Corresponding values are `"row"`, `"column"` and `"none"`

  `cluster_rows`

  :   boolean values determining if rows should be clustered or `hclust`
      object,

  `cluster_cols`

  :   boolean values determining if columns should be clustered or
      `hclust` object.

  `clustering_distance_rows`

  :   distance measure used in clustering rows. Possible values are
      `"correlation"` for Pearson correlation and all the distances
      supported by [`dist`](https://rdrr.io/r/stats/dist.html), such as
      `"euclidean"`, etc. If the value is none of the above it is
      assumed that a distance matrix is provided.

  `clustering_distance_cols`

  :   distance measure used in clustering columns. Possible values the
      same as for clustering_distance_rows.

  `clustering_method`

  :   clustering method used. Accepts the same values as
      [`hclust`](https://rdrr.io/r/stats/hclust.html).

  `clustering_callback`

  :   callback function to modify the clustering. Is called with two
      parameters: original `hclust` object and the matrix used for
      clustering. Must return a `hclust` object.

  `cutree_rows`

  :   number of clusters the rows are divided into, based on the
      hierarchical clustering (using cutree), if rows are not clustered,
      the argument is ignored

  `cutree_cols`

  :   similar to `cutree_rows`, but for columns

  `treeheight_row`

  :   the height of a tree for rows, if these are clustered. Default
      value 50 points.

  `treeheight_col`

  :   the height of a tree for columns, if these are clustered. Default
      value 50 points.

  `legend`

  :   logical to determine if legend should be drawn or not.

  `legend_breaks`

  :   vector of breakpoints for the legend.

  `legend_labels`

  :   vector of labels for the `legend_breaks`.

  `annotation_row`

  :   data frame that specifies the annotations shown on left side of
      the heatmap. Each row defines the features for a specific row. The
      rows in the data and in the annotation are matched using
      corresponding row names. Note that color schemes takes into
      account if variable is continuous or discrete.

  `annotation_col`

  :   similar to annotation_row, but for columns.

  `annotation`

  :   deprecated parameter that currently sets the annotation_col if it
      is missing

  `annotation_colors`

  :   list for specifying annotation_row and annotation_col track colors
      manually. It is possible to define the colors for only some of the
      features. Check examples for details.

  `annotation_legend`

  :   boolean value showing if the legend for annotation tracks should
      be drawn.

  `annotation_names_row`

  :   boolean value showing if the names for row annotation tracks
      should be drawn.

  `annotation_names_col`

  :   boolean value showing if the names for column annotation tracks
      should be drawn.

  `drop_levels`

  :   logical to determine if unused levels are also shown in the legend

  `show_rownames`

  :   boolean specifying if column names are be shown.

  `show_colnames`

  :   boolean specifying if column names are be shown.

  `main`

  :   the title of the plot

  `fontsize`

  :   base fontsize for the plot

  `fontsize_row`

  :   fontsize for rownames (Default: fontsize)

  `fontsize_col`

  :   fontsize for colnames (Default: fontsize)

  `angle_col`

  :   angle of the column labels, right now one can choose only from few
      predefined options (0, 45, 90, 270 and 315)

  `display_numbers`

  :   logical determining if the numeric values are also printed to the
      cells. If this is a matrix (with same dimensions as original
      matrix), the contents of the matrix are shown instead of original
      values.

  `number_format`

  :   format strings (C printf style) of the numbers shown in cells. For
      example "`%.2f`" shows 2 decimal places and "`%.1e`" shows
      exponential notation (see more in
      [`sprintf`](https://rdrr.io/r/base/sprintf.html)).

  `number_color`

  :   color of the text

  `fontsize_number`

  :   fontsize of the numbers displayed in cells

  `gaps_row`

  :   vector of row indices that show where to put gaps into heatmap.
      Used only if the rows are not clustered. See `cutree_row` to see
      how to introduce gaps to clustered rows.

  `gaps_col`

  :   similar to gaps_row, but for columns.

  `labels_row`

  :   custom labels for rows that are used instead of rownames.

  `labels_col`

  :   similar to labels_row, but for columns.

  `filename`

  :   file path where to save the picture. Filetype is decided by the
      extension in the path. Currently following formats are supported:
      png, pdf, tiff, bmp, jpeg. Even if the plot does not fit into the
      plotting window, the file size is calculated so that the plot
      would fit there, unless specified otherwise.

  `width`

  :   manual option for determining the output file width in inches.

  `height`

  :   manual option for determining the output file height in inches.

  `silent`

  :   do not draw the plot (useful when using the gtable output)

  `na_col`

  :   specify the color of the NA cell in the matrix.

## Examples

``` r
urt %>%
  aggregate_taxa(rank = "order") %>%
  tacoplot_prevalences(location, cutoff = .1,
  treeheight_row = 0, cutree_rows = 4,
  fontsize = 6, cellwidth = 15)

```
