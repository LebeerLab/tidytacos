# Add taxon prevalences to the taxon table

`add_prevalence()` calculates taxon prevalences (overall or per
condition) and adds it to the taxa table under the column name
"prevalence". Prevalence can be expressed as the number of samples where
a taxon occurs or the ratio of samples where a taxon occurs and the
total amount of samples.

## Usage

``` r
add_prevalence(ta, condition = NULL, relative = FALSE, fisher_test = FALSE)
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
