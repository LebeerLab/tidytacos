# tidytacos: Functions to manipulate and visualize amplicon count data.

tidytacos is an R package for the analysis of amplicon count data:
abundances of amplicon sequences (either clustered in OTUs or exact
variants) in samples. The package builds on the tidyverse created by
Hadley Wickham: the data are stored in tidy tables where each row is an
observation and each column a variable. In addition, the package
supplies a set of "verbs": functions that take a tidytacos object as
first argument and also return a tidytacos object. Not all functionality
is currently implemented in the form of verbs, but this will soon be
remediated.

## Author(s)

Stijn Wittouck <wittouck_stijn@hotmail.com>
