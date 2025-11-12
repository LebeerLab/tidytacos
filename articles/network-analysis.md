# Network Analysis

``` r
suppressMessages({
  library(tidyverse)
  library(tidytacos)
})

urt.nf <- urt %>% 
filter_samples(location == "NF") %>%
add_taxon_name() 
```

This is a brief guide on how to perform network analysis on microbiome
count data using a similar strategy as performed in the [Isala
study](https://www.nature.com/articles/s41564-023-01500-0) using
tidytacos.

The method clusters taxa based on their co-occurrence patterns across
the samples.

As per usual the test dataset used is a human microbiome samples from
the upper respiratory tract (URT), taken from [this
paper](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02372/full)
by De Boeck et al. It contains nose as well as nasopharynx samples. Most
samples were taken using a swab method, but a minority was taking with
the aspirate method.

### Step 1: Creating a network

The first step in the analysis is creating a sparcc network from our
taxonomic counts. For this the [SpiecEasi
package](https://github.com/zdk123/SpiecEasi) is used.

``` r
devtools::install_github('zdk123/SpiecEasi')
```

``` r
urt.net <- urt.nf %>%
  network(
    min_occurrence = 0.05, # the minimum percentage of samples that needs to have a taxa for it to be included in the network
    iter=20, # number of iterations of the outer loop of sparcc 
    inner_iter=10, # number of iterations of the inner loop of sparcc
    th=0.1, # correlations beneath this threshold are set to 0
    taxon_name=taxon_name)

# The network object is a list with two elements:
urt.net$Cor %>% head(2)# Correlation matrix
#>             [,1]        [,2]         [,3]        [,4]       [,5]       [,6]
#> [1,]  1.00000000 -0.02529089 -0.005801137  0.01516097 0.05986462 0.13709887
#> [2,] -0.02529089  1.00000000  0.030665238 -0.09163985 0.07760764 0.02597374
#>             [,7]        [,8]       [,9]       [,10]       [,11]       [,12]
#> [1,]  0.03071816 -0.08361307 0.10603959 -0.01524706 -0.09335080 -0.01788330
#> [2,] -0.05482966 -0.07955803 0.04296975 -0.01357689  0.08080531  0.05173014
#>            [,13]       [,14]       [,15]       [,16]       [,17]       [,18]
#> [1,] -0.01253595  0.03261072  0.01128624 -0.01741375 -0.01796809 0.308998423
#> [2,]  0.08546550 -0.03786864 -0.01991150 -0.03926070  0.02300305 0.006928278
#>            [,19]      [,20]        [,21]       [,22]      [,23]        [,24]
#> [1,] -0.06816382 0.04028286 -0.007349141 -0.01681414 0.04898641 -0.006463706
#> [2,] -0.08602481 0.08172153 -0.005358079  0.15721554 0.04186287  0.010783317
#>            [,25]       [,26]        [,27]       [,28]       [,29]       [,30]
#> [1,] -0.02932673 -0.04693648 -0.006408764 -0.08380906 -0.06849173 -0.09679495
#> [2,]  0.13554391 -0.10768510 -0.038097269  0.03864768  0.02522536  0.01495462
#>            [,31]      [,32]      [,33]       [,34]       [,35]       [,36]
#> [1,] 0.006154588 0.05529007  0.1787906 -0.01988641 -0.05646627  0.02435975
#> [2,] 0.124082370 0.12098392 -0.0169648 -0.07313655 -0.05252636 -0.01685275
#>            [,37]       [,38]       [,39]       [,40]       [,41]        [,42]
#> [1,]  0.05573063  0.27021927 -0.08233737  0.27717752  0.32687175 -0.020168237
#> [2,] -0.09939702 -0.04105809  0.10925131 -0.02920952 -0.04238792  0.004643286
#>            [,43]       [,44]       [,45]       [,46]       [,47]       [,48]
#> [1,] -0.13372090 -0.02110049 -0.11099158 -0.01791130 -0.05660303 -0.04285631
#> [2,] -0.03572358  0.04358548 -0.06357924 -0.02876454 -0.05476827  0.01373997
#>            [,49]      [,50]      [,51]      [,52]       [,53]        [,54]
#> [1,] -0.06429722 0.01189109 0.04833545 -0.0179242  0.01958685 -0.003913994
#> [2,] -0.01867380 0.09034714 0.01886988 -0.0573012 -0.01708747  0.026535551
#>             [,55]      [,56]        [,57]      [,58]      [,59]       [,60]
#> [1,] -0.011371461 0.09608058  0.008814129 0.08724145 0.02541135  0.05272446
#> [2,]  0.001016776 0.03030416 -0.003755814 0.09917971 0.06309768 -0.02141835
#>            [,61]       [,62]       [,63]       [,64]       [,65]       [,66]
#> [1,] 0.007878015 -0.07609307 -0.07515718 -0.08300833 -0.09037017  0.04006076
#> [2,] 0.072125421 -0.04605675 -0.09692066 -0.08869844  0.02120372 -0.06606351
#>            [,67]       [,68]       [,69]       [,70]       [,71]       [,72]
#> [1,] -0.07315554 -0.02608547 -0.02595041 -0.12857806 -0.02958168 -0.05685777
#> [2,]  0.11404219 -0.01836244  0.01069843 -0.05887034  0.04082737  0.08911415
#>            [,73]       [,74]      [,75]       [,76]      [,77]       [,78]
#> [1,] -0.02120009 0.015739045 0.03643760 0.006902477 0.01807948  0.14955167
#> [2,] -0.04430801 0.005226727 0.05359429 0.126778092 0.01243113 -0.09374607
#>           [,79]         [,80]      [,81]       [,82]        [,83]      [,84]
#> [1,] 0.32421126 -0.0203249517 0.06744001 -0.04116257 -0.007734448 0.09769376
#> [2,] 0.01576832 -0.0002167615 0.01841182  0.09279505  0.021557648 0.15094090
#>            [,85]      [,86]        [,87]       [,88]       [,89]        [,90]
#> [1,] -0.02065859 0.08824767 0.0002962855 -0.07892791 -0.01215166  0.023823958
#> [2,]  0.04044059 0.01583869 0.0379670151 -0.10003871 -0.01476276 -0.002405105
#>            [,91]       [,92]       [,93]       [,94]       [,95]       [,96]
#> [1,]  0.02127067 -0.04943964 -0.03150189  0.01665404 -0.02428586 -0.03945849
#> [2,] -0.04528266  0.10176461 -0.05882543 -0.04453790  0.05828724 -0.05052593
#>            [,97]       [,98]        [,99]      [,100]      [,101]      [,102]
#> [1,]  0.04364262 -0.02789056 -0.162288183 -0.06216937 -0.07105051 0.118253684
#> [2,] -0.01618929 -0.05109104 -0.007859769 -0.05034302 -0.14405057 0.001834282
#>           [,103]      [,104]     [,105]      [,106]      [,107]      [,108]
#> [1,] -0.04045566 -0.03360512 -0.0701020 -0.04941166 -0.01169133  0.02423498
#> [2,] -0.00955543 -0.07275669 -0.0178108 -0.04841629 -0.10043112 -0.07519695
#>           [,109]
#> [1,] -0.02405523
#> [2,] -0.04680215
# urt.net$Cov (covariance matrix)
# urt.net$names (names of the taxa on x and y-axes)
```

## Step 2: Filtering the network matrix

Next up we filter the output correlation matrix to remove weak
correlations. For this we use the `filter_network` function, which has
the `Matrix` package as a dependency. We also install the `igraph`
package, which can be used for the visualization of the network.

``` r
install.packages("Matrix")
install.packages("igraph")
```

We remove weak correlations (here set to an absolute value of \<=0.1)
and remove the diagonal. This is all handled with the `filter_network`
function.

``` r
urt.net.filt <- urt.net %>% filter_network(threshold = 0.1)
```

## Step 3: Clustering of the network

Now the taxa in the network are clustered using an unsupervised
clustering algorithm, the [Markov Cluster Algorithm
(MCL)](https://micans.org/mcl/). The `cluster_network` function has the
`MCL` package as a dependency. The algorithm cannot handle negative
interactions so any correlation \<0 will be set to 0.

``` r
install.packages("MCL")
```

``` r
urt.net.clust <- urt.net.filt %>% 
  cluster_network(min_n=5, visualize=T)
```

![](network-analysis_files/figure-html/network-clustering-1.png)

``` r

urt.net.clust %>% 
  group_by(cluster) %>% 
  summarize(n=n())
#> # A tibble: 5 × 2
#>   cluster     n
#>   <chr>   <int>
#> 1 c0         42
#> 2 c17         8
#> 3 c18         8
#> 4 c21        20
#> 5 c54        11
```

Using the new clustered taxa (“eigentaxa”) we can visualize the clusters
in the microbiome profile.

``` r
urt.nf$taxa <- urt.nf$taxa %>% 
  left_join(urt.net.clust, by=c("taxon_name"="taxon")) 
urt.nf <- urt.nf %>%
  add_taxon_name() %>% 
  mutate_taxa(cluster= ifelse(is.na(cluster), taxon_name, cluster)) %>% 
  set_rank_names(c("cluster")) %>% 
  aggregate_taxa(rank="cluster")
urt.nf  %>% tacoplot_stack()
```

![](network-analysis_files/figure-html/unnamed-chunk-2-1.png)

``` r

urt.net.clust %>% filter(cluster == unique(urt.net.clust$cluster)[[4]])
#> # A tibble: 20 × 2
#>    taxon            cluster
#>    <chr>            <chr>  
#>  1 Alloprevotella 3 c21    
#>  2 Campylobacter 2  c21    
#>  3 Fusobacterium 5  c21    
#>  4 Gemella          c21    
#>  5 Granulicatella 1 c21    
#>  6 Granulicatella 2 c21    
#>  7 Haemophilus 3    c21    
#>  8 Leptotrichia 1   c21    
#>  9 Megasphaera 1    c21    
#> 10 Neisseria 2      c21    
#> 11 Neisseria 4      c21    
#> 12 Prevotella_7 1   c21    
#> 13 Streptococcus 3  c21    
#> 14 Streptococcus 4  c21    
#> 15 Streptococcus 5  c21    
#> 16 Streptococcus 9  c21    
#> 17 Veillonella 1    c21    
#> 18 Veillonella 2    c21    
#> 19 Veillonella 4    c21    
#> 20 Veillonella 5    c21
```
