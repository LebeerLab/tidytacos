# Create extra variables in the count table

Create extra variables in the count table

## Usage

``` r
mutate_counts(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Mutate criteria for the counts table.

## Value

A tidytacos object.

## Examples

``` r
# add a column to the counts table
urt %>% mutate_counts(cl10 = log10(count))
#> $samples
#> # A tibble: 217 × 9
#>    run    condition participant location method plate passes_qc sample sample_id
#>    <chr>  <chr>     <chr>       <chr>    <chr>  <dbl> <lgl>     <chr>  <chr>    
#>  1 20161… CON       CON100      NF       S          3 TRUE      CON10… s1       
#>  2 20161… CON       CON100      N        S          3 TRUE      CON10… s2       
#>  3 20161… CON       CON10       NF       A          1 TRUE      CON10… s3       
#>  4 20161… CON       CON10       NF       S          1 TRUE      CON10… s4       
#>  5 20161… CON       CON11       NF       S          1 TRUE      CON11… s5       
#>  6 20161… CON       CON11       N        S          1 TRUE      CON11… s6       
#>  7 20161… CON       CON12       NF       S          1 TRUE      CON12… s7       
#>  8 20161… CON       CON13       NF       S          1 TRUE      CON13… s8       
#>  9 20161… CON       CON13       N        S          1 TRUE      CON13… s9       
#> 10 20161… CON       CON14       NF       S          1 TRUE      CON14… s10      
#> # ℹ 207 more rows
#> 
#> $taxa
#> # A tibble: 1,957 × 10
#>    taxon       kingdom phylum class order family genus species sequence taxon_id
#>    <chr>       <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>    <chr>   
#>  1 TACAGAGGGT… Bacter… Prote… Gamm… Pseu… Morax… Mora… catarr… TACAGAG… t1      
#>  2 TACGTAGGTG… Bacter… Firmi… Baci… Baci… Staph… Stap… aureus… TACGTAG… t2      
#>  3 TACAGAGGGT… Bacter… Prote… Gamm… Pseu… Morax… Mora… porci   TACAGAG… t3      
#>  4 TACAGAGGGT… Bacter… Prote… Gamm… Pseu… Morax… Mora… bovis/… TACAGAG… t4      
#>  5 TACGTAGGGT… Bacter… Prote… Beta… Neis… Neiss… Neis… mening… TACGTAG… t5      
#>  6 TACGTATGTC… Bacter… Fusob… Fuso… Fuso… Fusob… Fuso… canife… TACGTAT… t6      
#>  7 TACGTAGGGT… Bacter… Actin… Acti… Cory… Coryn… Cory… accole… TACGTAG… t7      
#>  8 TACGGAGGGT… Bacter… Prote… Gamm… Past… Paste… Haem… haemol… TACGGAG… t8      
#>  9 TACGTATGTC… Bacter… Fusob… Fuso… Fuso… Fusob… Fuso… nuclea… TACGTAT… t9      
#> 10 TACGTAGGGT… Bacter… Prote… Beta… Neis… Neiss… Neis… lactam… TACGTAG… t10     
#> # ℹ 1,947 more rows
#> 
#> $counts
#> # A tibble: 7,693 × 4
#>    count sample_id taxon_id  cl10
#>    <dbl> <chr>     <chr>    <dbl>
#>  1 13473 s1        t7        4.13
#>  2 10322 s2        t7        4.01
#>  3    86 s4        t7        1.93
#>  4 19770 s5        t7        4.30
#>  5 28209 s6        t7        4.45
#>  6   416 s7        t7        2.62
#>  7  9508 s8        t7        3.98
#>  8 31623 s9        t7        4.50
#>  9   130 s10       t7        2.11
#> 10   922 s11       t7        2.96
#> # ℹ 7,683 more rows
#> 
#> $rank_names
#> [1] "kingdom" "phylum"  "class"   "order"   "family"  "genus"  
#> 
#> attr(,"class")
#> [1] "tidytacos"
```
