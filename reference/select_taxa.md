# Retain or remove a set of taxon variables

Retain or remove a set of taxon variables

## Usage

``` r
select_taxa(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Selection criteria for the taxa table.

## Value

A tidytacos object.

## Examples

``` r
# drop the sequence column
urt %>% select_taxa(-sequence)
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
#> # A tibble: 1,957 × 9
#>    taxon                kingdom phylum class order family genus species taxon_id
#>    <chr>                <chr>   <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr>   
#>  1 TACAGAGGGTGCAAGCGTT… Bacter… Prote… Gamm… Pseu… Morax… Mora… catarr… t1      
#>  2 TACGTAGGTGGCAAGCGTT… Bacter… Firmi… Baci… Baci… Staph… Stap… aureus… t2      
#>  3 TACAGAGGGTGCAAGCGTT… Bacter… Prote… Gamm… Pseu… Morax… Mora… porci   t3      
#>  4 TACAGAGGGTGCAAGCGTT… Bacter… Prote… Gamm… Pseu… Morax… Mora… bovis/… t4      
#>  5 TACGTAGGGTGCGAGCGTT… Bacter… Prote… Beta… Neis… Neiss… Neis… mening… t5      
#>  6 TACGTATGTCACAAGCGTT… Bacter… Fusob… Fuso… Fuso… Fusob… Fuso… canife… t6      
#>  7 TACGTAGGGTGCGAGCGTT… Bacter… Actin… Acti… Cory… Coryn… Cory… accole… t7      
#>  8 TACGGAGGGTGCGAGCGTT… Bacter… Prote… Gamm… Past… Paste… Haem… haemol… t8      
#>  9 TACGTATGTCACGAGCGTT… Bacter… Fusob… Fuso… Fuso… Fusob… Fuso… nuclea… t9      
#> 10 TACGTAGGGTGCGAGCGTT… Bacter… Prote… Beta… Neis… Neiss… Neis… lactam… t10     
#> # ℹ 1,947 more rows
#> 
#> $counts
#> # A tibble: 7,693 × 3
#>    count sample_id taxon_id
#>    <dbl> <chr>     <chr>   
#>  1 13473 s1        t7      
#>  2 10322 s2        t7      
#>  3    86 s4        t7      
#>  4 19770 s5        t7      
#>  5 28209 s6        t7      
#>  6   416 s7        t7      
#>  7  9508 s8        t7      
#>  8 31623 s9        t7      
#>  9   130 s10       t7      
#> 10   922 s11       t7      
#> # ℹ 7,683 more rows
#> 
#> $rank_names
#> [1] "kingdom" "phylum"  "class"   "order"   "family"  "genus"  
#> 
#> attr(,"class")
#> [1] "tidytacos"

# keep only the taxon_id and genus columns
urt %>% select_taxa(taxon_id, genus)
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
#> # A tibble: 1,957 × 2
#>    taxon_id genus            
#>    <chr>    <chr>            
#>  1 t1       Moraxella        
#>  2 t2       Staphylococcus   
#>  3 t3       Moraxella        
#>  4 t4       Moraxella        
#>  5 t5       Neisseria        
#>  6 t6       Fusobacterium    
#>  7 t7       Corynebacterium_1
#>  8 t8       Haemophilus      
#>  9 t9       Fusobacterium    
#> 10 t10      Neisseria        
#> # ℹ 1,947 more rows
#> 
#> $counts
#> # A tibble: 7,693 × 3
#>    count sample_id taxon_id
#>    <dbl> <chr>     <chr>   
#>  1 13473 s1        t7      
#>  2 10322 s2        t7      
#>  3    86 s4        t7      
#>  4 19770 s5        t7      
#>  5 28209 s6        t7      
#>  6   416 s7        t7      
#>  7  9508 s8        t7      
#>  8 31623 s9        t7      
#>  9   130 s10       t7      
#> 10   922 s11       t7      
#> # ℹ 7,683 more rows
#> 
#> $rank_names
#> [1] "genus"
#> 
#> attr(,"class")
#> [1] "tidytacos"
```
