% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{add_eigentaxa}
\alias{add_eigentaxa}
\title{Calculates eigentaxa values based on SparCC - MCL generated clusters per sample. 
It is advised to run \code{\link{cluster_taxa}} on the tidyamplicons object
before running this function to add the clusters if you want to stray from any default parameters.}
\usage{
add_eigentaxa(ta, taxon_name = taxon, sample_name = sample)
}
\arguments{
\item{ta}{a tidytacos object.}
}
\description{
Calculates eigentaxa values based on SparCC - MCL generated clusters per sample. 
It is advised to run \code{\link{cluster_taxa}} on the tidyamplicons object
before running this function to add the clusters if you want to stray from any default parameters.
}
