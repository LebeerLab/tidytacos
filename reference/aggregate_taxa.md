# Aggregate taxa on a given taxonomic rank

There are two ways to call this function:

## Usage

``` r
aggregate_taxa(ta, rank = NULL)
```

## Arguments

- ta:

  A tidytacos object.

- rank:

  An optional rank to aggregate on.

## Value

A tidytacos object.

## Details

- If the rank you are interested in is in the standard list, just supply
  it as an argument.

- If not, delete all taxon variables except taxon_id and the ranks you
  are still interested in prior to calling this function.

## Examples

``` r
urt %>% aggregate_taxa(rank = "class")
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
#> # A tibble: 65 × 6
#>    kingdom  phylum         class                 taxon_id taxon_name       taxon
#>    <chr>    <chr>          <chr>                 <chr>    <chr>            <chr>
#>  1 Bacteria Firmicutes     Bacilli               t2       Bacilli          Baci…
#>  2 Bacteria Actinobacteria Actinobacteria        t5       Actinobacteria   Acti…
#>  3 Bacteria Proteobacteria Gammaproteobacteria   t1       Gammaproteobact… Gamm…
#>  4 Bacteria Firmicutes     Clostridia            t6       Clostridia       Clos…
#>  5 Bacteria Proteobacteria Betaproteobacteria    t3       Betaproteobacte… Beta…
#>  6 Bacteria Fusobacteria   Fusobacteriia         t4       Fusobacteriia    Fuso…
#>  7 Bacteria Firmicutes     Negativicutes         t7       Negativicutes    Nega…
#>  8 Bacteria Bacteroidetes  Bacteroidia           t8       Bacteroidia      Bact…
#>  9 Bacteria Proteobacteria Alphaproteobacteria   t9       Alphaproteobact… Alph…
#> 10 Bacteria Proteobacteria Epsilonproteobacteria t11      Epsilonproteoba… Epsi…
#> # ℹ 55 more rows
#> 
#> $counts
#> # A tibble: 2,064 × 3
#>    taxon_id sample_id count
#>    <chr>    <chr>     <dbl>
#>  1 t1       s1           85
#>  2 t1       s10       15606
#>  3 t1       s100       1103
#>  4 t1       s102       2777
#>  5 t1       s103        282
#>  6 t1       s104        587
#>  7 t1       s105        376
#>  8 t1       s106        118
#>  9 t1       s107       1561
#> 10 t1       s108       1628
#> # ℹ 2,054 more rows
#> 
#> $rank_names
#> [1] "kingdom" "phylum"  "class"  
#> 
#> attr(,"class")
#> [1] "tidytacos"
```
