# Trim all sequences

`trim_asvs()` trims sequence ends of the sequence supplied in the taxa
table. This function assumes that the sequence variable in the taxon
table is called "sequence".

## Usage

``` r
trim_asvs(ta, start, end)
```

## Arguments

- ta:

  A tidytacos object.

- start:

  Index of where to start trimming.

- end:

  Index of where to stop trimming.

## Value

A tidytacos object.

## Examples

``` r
# keep only the first 200 nucleotides of the sequences
urt %>% trim_asvs(0, 200)
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
#> # A tibble: 1,854 × 8
#>    taxon_id genus                   family   order class phylum kingdom sequence
#>    <chr>    <chr>                   <chr>    <chr> <chr> <chr>  <chr>   <chr>   
#>  1 t1       Moraxella               Moraxel… Pseu… Gamm… Prote… Bacter… TACAGAG…
#>  2 t100     Staphylococcus          Staphyl… Baci… Baci… Firmi… Bacter… TACGTAG…
#>  3 t1000    Pseudarcicella          Cytopha… Cyto… Cyto… Bacte… Bacter… TACGGAG…
#>  4 t1002    Leptotrichia            Leptotr… Fuso… Fuso… Fusob… Bacter… TACGTAT…
#>  5 t1003    hgcI_clade              Sporich… Fran… Acti… Actin… Bacter… TACATAG…
#>  6 t1004    Prevotella              Prevote… Bact… Bact… Bacte… Bacter… TACGGAA…
#>  7 t1005    Peptoniphilus           Family_… Clos… Clos… Firmi… Bacter… TACGTAG…
#>  8 t1006    Hymenobacter            Cytopha… Cyto… Cyto… Bacte… Bacter… TACGGAG…
#>  9 t1008    Dialister               Veillon… Sele… Nega… Firmi… Bacter… TACGTAG…
#> 10 t1009    Lachnospiraceae_UCG-005 Lachnos… Clos… Clos… Firmi… Bacter… TACGTAT…
#> # ℹ 1,844 more rows
#> 
#> $counts
#> # A tibble: 7,693 × 3
#>    sample_id taxon_id count
#>    <chr>     <chr>    <dbl>
#>  1 s1        t108        60
#>  2 s1        t11       2703
#>  3 s1        t14       5661
#>  4 s1        t146       437
#>  5 s1        t1511       21
#>  6 s1        t1792        8
#>  7 s1        t2        4153
#>  8 s1        t203         8
#>  9 s1        t23         83
#> 10 s1        t30       2646
#> # ℹ 7,683 more rows
#> 
#> $rank_names
#> [1] "kingdom" "phylum"  "class"   "order"   "family"  "genus"  
#> 
#> attr(,"class")
#> [1] "tidytacos"
```
