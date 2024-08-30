tidytacos <img src="man/figures/logo.png" align="right" width="200"/>
======================
[![R-CMD-check](https://github.com/LebeerLab/tidytacos/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LebeerLab/tidytacos/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/LebeerLab/tidytacos/graph/badge.svg?token=532YS16DXU)](https://codecov.io/gh/LebeerLab/tidytacos)
<br><br><br><br><br>
## Overview

Tidytacos (tidy TAxonomic COmpositionS) is an R package for the exploration of microbial community data. Such community data consists of read counts generated by amplicon sequencing (e.g. a region of the 16S rRNA gene) or metagenome (shotgun) sequencing. Each read count represents a number of sequencing reads identified for some taxon (an ASV, OTU, species, or higher-level taxon) in a sample. 

Tidytacos builds on the [tidyverse](https://www.tidyverse.org/) created by [Hadley Wickham](http://hadley.nz/): the data are stored in tidy tables where each row is an observation and each column a variable. In addition, the package supplies a set of "verbs": functions that take a tidytacos object as first argument and also return a tidytacos object. This makes it easy to construct "pipe chains" of code that represent series of operations performed on the tidytacos object. 

## Prerequisites 

Tidytacos is an R package. You can find instructions to download and install R [here](https://cran.r-project.org/).

Tidytacos relies on the tidyverse R package (or, more accurately, set of R packages). You can install the tidyverse by running the following R code: 

```R
install.packages("tidyverse")
```

Finally, RStudio is a nice IDE to work with R code (as well as code in other scripting languages). It has a lot more features than what the default R IDE allows: beyond creating and saving scripts, it also shows your figures, allows you to navigate files, allows you to inspect tables etc. You can download RStudio [here](https://posit.co/downloads/). 

## Installation

Run the following R code to install the latest version of tidytacos: 

```R
install.packages("devtools")
devtools::install_github("LebeerLab/tidytacos")
```

## Getting started

A tidytacos object is read and stored as three sparse tables (counts-, taxa- and samples.csv). 
To read in existing data from a folder, [for example one called ‘leaf’ in the ‘data-raw/tidytacos’ folder](https://github.com/LebeerLab/tidytacos/tree/dev/data-raw/tidytacos/leaf) you would run:

```R
taco <- read_tidytacos("data-raw/tidytacos/leaf")
```
If you have data in the form of a phyloseq object you could convert it using:

```R
taco <- from_phyloseq(phylo_obj)
```
If your ASVs are counted and annotated using [dada2](https://benjjneb.github.io/dada2/), you can use the following function to convert the results to a tidytacos object:
```R
taco <- from_dada(seqtab.nochim, taxa)
```
Where seqtab.nochim and taxa refer to the R objects [as calculated in the dada2 tutorial](https://benjjneb.github.io/dada2/tutorial.html)

Note that the taxa table includes all taxonomic levels from kingdom to species, in addition to the optional sequence variable (which refers to the nucleotide sequence) and the taxon_id variable which links the taxa table with the counts table in the tidytacos object.


## Documentation

[A documentation page (help page)](https://lebeerlab.github.io/tidytacos/reference/index.html) is available for all functions in the browser or in R. You can view it in R by running e.g. `?filter_samples`. Some useful tutorials can be found on the [wiki](https://github.com/LebeerLab/tidytacos/wiki). 

## Need support?

Post on [GitHub issues](https://github.com/LebeerLab/tidytacos/issues) if you have questions, requests, or if you run into an issue.

## Feel like contributing?

Please read the [GitHub Developer Guide](https://github.com/LebeerLab/tidytacos/wiki/Developer-Guide). Fork the dev branch, make your changes and make a pull request. Your suggestions will be reviewed and if approved, will be implemented in the next release.

