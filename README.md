tidytacos <img src="man/figures/logo.png" align="right" width="200"/>
======================
[![R-CMD-check](https://github.com/LebeerLab/tidytacos/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LebeerLab/tidytacos/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/LebeerLab/tidytacos/branch/feature/unit-testing/graph/badge.svg?token=WK0QN34DJB)](https://codecov.io/github/LebeerLab/tidytacos)

## Overview

tidytacos (tidy TAxonomic COmpositionS) is an R package for the exploration of microbial community data. Such community data consists of read counts generated by amplicon sequencing (e.g. a region of the 16S rRNA gene) or metagenome (shotgun) sequencing. Each read count represents a number of sequencing reads identified for some taxon (an OTU, ASV, species, or higher-level taxon) in a sample. 

tidytacos builds on the [tidyverse](https://www.tidyverse.org/) created by [Hadley Wickham](http://hadley.nz/): the data are stored in tidy tables where each row is an observation and each column a variable. In addition, the package supplies a set of "verbs": functions that take a tidytacos object as first argument and also return a tidytacos object.

## Prerequisites 

tidytacos is an R package. You can find instructions to download and install R [here](https://cran.r-project.org/).

tidytacos relies on the tidyverse R package (or, more accurately, set of R packages). You can install the tidyverse by running the following R code: 

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

## Documentation

A documentation page (help page) is available for all functions. You can view it by running e.g. `?filter_samples`. For most functions it is still very basic; this will be improved in the future. Some useful tutorials can be found on the [wiki](https://github.com/LebeerLab/tidytacos/wiki). 
