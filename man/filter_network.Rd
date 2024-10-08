% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{filter_network}
\alias{filter_network}
\title{Filters the output of \code{\link[=network]{network()}} to a minimal threshold
and transforms to matrix for downstream clustering or heatplot visualization.}
\usage{
filter_network(network, threshold = 0.1, fdr = NULL, ...)
}
\arguments{
\item{network}{network object generated by \code{\link[=network]{network()}}.}

\item{threshold}{absolute value of correlations below this threshold are}

\item{fdr}{the threshold for false discovery rate
(if pvalues are calculated for the network) filtered out.}

\item{...}{
  Arguments passed on to \code{\link[stats:p.adjust]{stats::p.adjust}}
  \describe{
    \item{\code{p}}{numeric vector of p-values (possibly with
    \code{\link{NA}}s).  Any other \R object is coerced by
    \code{\link{as.numeric}}.}
    \item{\code{method}}{correction method, a \code{\link{character}} string.
    Can be abbreviated.}
    \item{\code{n}}{number of comparisons, must be at least \code{length(p)};
    only set this (to non-default) when you know what you are doing!}
  }}
}
\description{
Filters the output of \code{\link[=network]{network()}} to a minimal threshold
and transforms to matrix for downstream clustering or heatplot visualization.
}
