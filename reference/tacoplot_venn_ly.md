# Return an interactive venn diagram of overlapping taxon_ids between conditions

Return an interactive venn diagram of overlapping taxon_ids between
conditions

## Usage

``` r
tacoplot_venn_ly(ta, condition, ...)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  The name of a variable in the samples table that contains a
  categorical value.

- ...:

  Arguments passed on to
  [`ggVennDiagram::ggVennDiagram`](https://gaospecial.github.io/ggVennDiagram/reference/ggVennDiagram.html)

  `x`

  :   list of items

  `category.names`

  :   default is names(x)

  `show_intersect`

  :   if TRUE the text can be visualized by \`plotly\`

  `set_color`

  :   color of set labels ("black")

  `set_size`

  :   size of set labels (NA)

  `label`

  :   format of region labels, select one from
      c("count","percent","both","none")

  `label_alpha`

  :   set 0 to remove the background of region labels

  `label_geom`

  :   layer of region labels, choose from c("label", "text")

  `label_color`

  :   color of region labels ("black")

  `label_size`

  :   size of region labels (NA)

  `label_percent_digit`

  :   number of digits when formatting percent label (0)

  `label_txtWidth`

  :   width of text used in showing intersect members, will be ignored
      unless show_intersection is TRUE (40)

  `edge_lty`

  :   line type of set edges ("solid")

  `edge_size`

  :   line width of set edges (1)

  `force_upset`

  :   if TRUE, will always produce Upset plot no matter how many sets
      have (FALSE)

  `nintersects`

  :   number of intersects. If NULL, all intersections will show.

  `order.intersect.by`

  :   'size', 'name', or "none"

  `order.set.by`

  :   'size', 'name', or "none"

  `relative_height`

  :   the relative height of top panel in upset plot

  `relative_width`

  :   the relative width of left panel in upset plot
