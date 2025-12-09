# Add taxon prevalences to the taxon table

`add_prevalence()` calculates taxon prevalences (overall or per
condition) and adds it to the taxa table under the column name
"prevalence". Prevalence can be expressed as the number of samples where
a taxon occurs or the ratio of samples where a taxon occurs and the
total amount of samples.

## Usage

``` r
add_prevalence(
  ta,
  condition = NULL,
  relative = FALSE,
  fisher_test = FALSE,
  ...
)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  A categorical variable (string).

- relative:

  Whether to use relative occurrences.

- fisher_test:

  Whether to perform a fisher test and add the p-values of the test to
  the taxa table.

- ...:

  Arguments passed on to
  [`stats::fisher.test`](https://rdrr.io/r/stats/fisher.test.html)

  `x`

  :   either a two-dimensional contingency table in matrix form, or a
      factor object.

  `y`

  :   a factor object; ignored if `x` is a matrix.

  `workspace`

  :   an integer specifying the size of the workspace used in the
      network algorithm. In units of 4 bytes. Only used for
      non-simulated p-values larger than \\2 \times 2\\ tables. Since R
      version 3.5.0, this also increases the internal stack size which
      allows larger problems to be solved, however sometimes needing
      hours. In such cases, `simulate.p.values=TRUE` may be more
      reasonable.

  `hybrid`

  :   a logical. Only used for larger than \\2 \times 2\\ tables, in
      which cases it indicates whether the exact probabilities (default)
      or a hybrid approximation thereof should be computed.

  `hybridPars`

  :   a numeric vector of length 3, by default describing “Cochran's
      conditions” for the validity of the chi-squared approximation, see
      ‘Details’.

  `control`

  :   a list with named components for low level algorithm control. At
      present the only one used is `"mult"`, a positive integer \\\ge
      2\\ with default 30 used only for larger than \\2 \times 2\\
      tables. This says how many times as much space should be allocated
      to paths as to keys: see file `fexact.c` in the sources of this
      package.

  `or`

  :   the hypothesized odds ratio. Only used in the \\2 \times 2\\ case.

  `alternative`

  :   indicates the alternative hypothesis and must be one of
      `"two.sided"`, `"greater"` or `"less"`. You can specify just the
      initial letter. Only used in the \\2 \times 2\\ case.

  `conf.int`

  :   logical indicating if a confidence interval for the odds ratio in
      a \\2 \times 2\\ table should be computed (and returned).

  `conf.level`

  :   confidence level for the returned confidence interval. Only used
      in the \\2 \times 2\\ case and if `conf.int = TRUE`.

  `simulate.p.value`

  :   a logical indicating whether to compute p-values by Monte Carlo
      simulation, in larger than \\2 \times 2\\ tables.

  `B`

  :   an integer specifying the number of replicates used in the Monte
      Carlo test.

## Value

A tidytacos object.

## Details

If 'condition' is specified, the prevalences will be calculated
separately for each group defined by the condition variable. This
variable should be present in the sample table.

If `condition` is specified, differential prevalence testing can be
performed by setting the `fisher_test` argument. Options are F (default)
or T. When set to T, significance of differential prevalence will be
added to the taxa table under column name `fisher_p`.

Condition should be a categorical variable present in the samples table.
Supply condition as a string.

## See also

Other taxa-modifiers:
[`add_mean_rel_abundance()`](https://lebeerlab.github.io/tidytacos/reference/add_mean_rel_abundance.md),
[`add_taxon_name()`](https://lebeerlab.github.io/tidytacos/reference/add_taxon_name.md),
[`add_taxon_name_color()`](https://lebeerlab.github.io/tidytacos/reference/add_taxon_name_color.md)

## Examples

``` r
# add prevalences of all taxa
urtp <- urt %>% add_prevalence()
urtp$taxa %>% dplyr::select(taxon_id, prevalence)
#> # A tibble: 1,957 × 2
#>    taxon_id prevalence
#>    <chr>         <int>
#>  1 t1               86
#>  2 t2              207
#>  3 t3                7
#>  4 t4                5
#>  5 t5               16
#>  6 t6               22
#>  7 t7              155
#>  8 t8               25
#>  9 t9               24
#> 10 t10               1
#> # ℹ 1,947 more rows

# add prevalences and fisher test for location
urtpf <- urt %>%
  add_prevalence(condition="location", fisher_test=TRUE, relative=TRUE)
urtpf$taxa %>%
  dplyr::select(taxon_id, prevalence_in_N, prevalence_in_NF, fisher_p)
#> # A tibble: 1,957 × 4
#>    taxon_id prevalence_in_N prevalence_in_NF fisher_p
#>    <chr>              <dbl>            <dbl>    <dbl>
#>  1 t1                0.326           0.443    0.0885 
#>  2 t2                0.988           0.931    0.0437 
#>  3 t3                0.0233          0.0382   0.706  
#>  4 t4                0.0233          0.0229   1      
#>  5 t5                0.0581          0.0840   0.599  
#>  6 t6                0.128           0.0840   0.359  
#>  7 t7                0.814           0.649    0.0120 
#>  8 t8                0.0465          0.160    0.00933
#>  9 t9                0.0349          0.160    0.00354
#> 10 t10               0               0.00763  1      
#> # ℹ 1,947 more rows
```
