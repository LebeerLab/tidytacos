# Return a scree plot to visualize the eigenvalues of the PCA.

Return a scree plot to visualize the eigenvalues of the PCA.

## Usage

``` r
tacoplot_scree(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Arguments passed on to
  [`factoextra::fviz_eig`](https://rdrr.io/pkg/factoextra/man/eigenvalue.html)

  `X`

  :   an object of class PCA, CA, MCA, FAMD, MFA and HMFA
      \[FactoMineR\]; prcomp and princomp \[stats\]; dudi, pca, coa and
      acm \[ade4\]; ca and mjca \[ca package\].

  `choice`

  :   a text specifying the data to be plotted. Allowed values are
      "variance" or "eigenvalue".

  `geom`

  :   a text specifying the geometry to be used for the graph. Allowed
      values are "bar" for barplot, "line" for lineplot or c("bar",
      "line") to use both types.

  `barfill`

  :   fill color for bar plot.

  `barcolor`

  :   outline color for bar plot.

  `linecolor`

  :   color for line plot (when geom contains "line").

  `ncp`

  :   a numeric value specifying the number of dimensions to be shown.

  `addlabels`

  :   logical value. If TRUE, labels are added at the top of bars or
      points showing the information retained by each dimension.

  `hjust`

  :   horizontal adjustment of the labels.

  `main,xlab,ylab`

  :   plot main and axis titles.

  `ggtheme`

  :   function, ggplot2 theme name. Default value is theme_pubr().
      Allowed values include ggplot2 official themes: theme_gray(),
      theme_bw(), theme_minimal(), theme_classic(), theme_void(), ....

## Examples

``` r
urt %>% tacoplot_scree()
#> Warning: Ignoring empty aesthetic: `width`.
```
