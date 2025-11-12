# Retain or remove a set of count variables

Retain or remove a set of count variables

## Usage

``` r
select_counts(ta, ...)
```

## Arguments

- ta:

  A tidytacos object.

- ...:

  Selection criteria for the counts table.

## Value

A tidytacos object.

## Examples

``` r
# add a column to the counts table
leaf_ab <- leaf %>% add_rel_abundance()
# remove that column again
leaf_ab %>% select_counts(-rel_abundance)
#> $samples
#> # A tibble: 33 × 13
#>    sample description     Plant Spraydate Plot  Samplingdate Timepoint   Day
#>    <chr>  <chr>           <chr> <chr>     <chr> <date>           <dbl> <dbl>
#>  1 S224   BLANK-0-KIT1    Blank NA        KIT1  NA                  NA    NA
#>  2 S84    BLANK-1-S-KIT1  Blank NA        KIT1  NA                  NA    NA
#>  3 S220   BLANK-10-S-KIT3 Blank NA        KIT3  NA                  NA    NA
#>  4 S40    BLANK-2-KIT1    Blank NA        KIT1  NA                  NA    NA
#>  5 S89    BLANK-3-S-KIT1  Blank NA        KIT1  NA                  NA    NA
#>  6 S22    BLANK-4-S-KIT1  Blank NA        KIT1  NA                  NA    NA
#>  7 S35    BLANK-5-S-KIT2  Blank NA        KIT2  NA                  NA    NA
#>  8 S98    BLANK-6-KIT2    Blank NA        KIT2  NA                  NA    NA
#>  9 S124   BLANK-7-S-KIT2  Blank NA        KIT2  NA                  NA    NA
#> 10 S189   BLANK-7b-S-KIT2 Blank NA        KIT2  NA                  NA    NA
#> # ℹ 23 more rows
#> # ℹ 5 more variables: Leafweight <dbl>, Square <chr>, Treatment <chr>,
#> #   added_spike_copies <dbl>, sample_id <chr>
#> 
#> $taxa
#> # A tibble: 6,054 × 9
#>    kingdom  phylum         class       order family genus species taxon taxon_id
#>    <chr>    <chr>          <chr>       <chr> <chr>  <chr> <chr>   <chr> <chr>   
#>  1 Eukarya  Streptophyta   eudicotyle… core… Convo… Cusc… NA      TACA… t1      
#>  2 Eukarya  Streptophyta   Liliopsida  Poal… Poace… Agro… NA      GACA… t2      
#>  3 Eukarya  Streptophyta   Liliopsida  NA    NA     NA    NA      GACG… t3      
#>  4 Eukarya  Streptophyta   eudicotyle… core… Brass… Bras… napus   AACG… t4      
#>  5 Eukarya  Streptophyta   eudicotyle… core… Convo… Cusc… NA      TACA… t5      
#>  6 Bacteria Proteobacteria Gammaprote… Ente… Enter… Citr… NA      TACG… t6      
#>  7 Eukarya  Streptophyta   eudicotyle… core… Convo… Cusc… NA      GACA… t7      
#>  8 Bacteria Proteobacteria Gammaprote… Ente… Enter… Buch… NA      TACG… t8      
#>  9 Bacteria Proteobacteria Betaproteo… Burk… Oxalo… Mass… aurea/… TACG… t9      
#> 10 Eukarya  Oomycota       Oomycetes   Pyth… Pythi… Pyth… NA      TACG… t10     
#> # ℹ 6,044 more rows
#> 
#> $counts
#> # A tibble: 9,191 × 3
#>    count sample_id taxon_id
#>    <dbl> <chr>     <chr>   
#>  1    15 s1        t1      
#>  2    34 s2        t1      
#>  3    14 s3        t1      
#>  4    20 s4        t1      
#>  5    12 s5        t1      
#>  6    11 s6        t1      
#>  7    17 s7        t1      
#>  8    51 s8        t1      
#>  9    29 s9        t1      
#> 10    33 s10       t1      
#> # ℹ 9,181 more rows
#> 
#> $rank_names
#> [1] "kingdom" "phylum"  "order"   "class"   "family"  "genus"  
#> 
#> attr(,"class")
#> [1] "tidytacos"
```
