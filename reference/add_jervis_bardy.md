# Apply the taxon QC method of Jervis-Bardy

`add_jervis_bardy()` calculates the spearman correlation between
relative abundance and sample DNA concentration, for each taxon and adds
the correlation metric and p-value to the taxa table under the column
names "jb_cor" and "jb_p", respectively. If taxa show a distribution
that is negatively correlated with DNA concentration, it indicates their
potential as contaminants.

## Usage

``` r
add_jervis_bardy(ta, dna_conc, sample_condition = TRUE, min_pres = 3)
```

## Arguments

- ta:

  A tidytacos object.

- dna_conc:

  A variable in the samples table that contains dna concetrations
  (unquoted).

- sample_condition:

  An optional extra condition that samples must pass before
  calculations.

- min_pres:

  The minimum number of samples a taxon has to be present in for its
  correlation to be calculated.

## Value

A tidytacos object with the Jervis-Bardy metrics added to the taxa
table.

## Details

See: J. Jervis-Bardy et al., “Deriving accurate microbiota profiles from
human samples with low bacterial content through post-sequencing
processing of Illumina MiSeq data,” Microbiome, vol. 3, no. 1, Art. no.
1, 2015, doi: 10.1186/s40168-015-0083-8.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following object is masked from ‘package:tidytacos’:
#> 
#>     everything
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
# filter out blank samples
plants <- leaf %>%
  filter_samples(Plant != "Blank")
# assume Leafweight is a proxy for DNA concentration of the sample
plants_jb <- plants %>%
  add_jervis_bardy(dna_conc = Leafweight)

# we can do this in one step!
plants_jb <- leaf %>%
  add_jervis_bardy(
    dna_conc = Leafweight,
    sample_condition = Plant != "Blank"
)

# show the negative correlations
plants_jb$taxa %>%
  select(taxon_id, starts_with("jb_")) %>%
  filter(jb_cor < 0) %>%
  arrange(jb_p)
#> # A tibble: 409 × 3
#>    taxon_id jb_cor    jb_p
#>    <chr>     <dbl>   <dbl>
#>  1 t256     -0.964 0.00139
#>  2 t133     -0.706 0.00664
#>  3 t21      -0.582 0.0100 
#>  4 t212     -0.886 0.0167 
#>  5 t471     -0.886 0.0167 
#>  6 t30      -0.521 0.0295 
#>  7 t195     -0.533 0.0321 
#>  8 t111     -1     0.0417 
#>  9 t199     -1     0.0417 
#> 10 t279     -1     0.0417 
#> # ℹ 399 more rows
```
