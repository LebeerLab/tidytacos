# Perform compositional differential abundance analysis

`add_codifab()` performs a differential abundance test for all pairwise
ratios between taxa. Taxa that have a relatively high number of
significantly different ratios, can be considered more abundant in one
condition versus the other. The
[`tacoplot_codifab()`](https://lebeerlab.github.io/tidytacos/reference/tacoplot_codifab.md)
function allows better interpretation of these results.

## Usage

``` r
add_codifab(ta, condition, conditions = NULL, max_taxa = 30)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  A binary variable in the sample table (unquoted).

- conditions:

  A character vector with exactly two categories of the condition
  variable.

- max_taxa:

  The maximum number of taxa to use.

## Value

A tidytacos object with an extra table taxon_pairs

## Details

A table called taxon_pairs will be added to the tidytacos object, with
for each pair of a taxon and a reference taxon, the differential
abundance of the taxon between the two conditions (with respect to the
reference taxon). The test that is performed is a Wilcoxon rank sum test
and the test statistic that is reported is the two-sample Hodges–Lehmann
estimator (the median of all pairwise differences between the samples).

It is possible to supply the conditions to compare through the
conditions argument. Other conditions than the two supplied will be
removed from the data.

This method is based on the principle introduced by Aitchison in "The
statistical analysis of compositional data." Journal of the Royal
Statistical Society: Series B (Methodological) 44.2 (1982): 139-16

## See also

Other codifab-functions:
[`tacoplot_codifab()`](https://lebeerlab.github.io/tidytacos/reference/tacoplot_codifab.md)
