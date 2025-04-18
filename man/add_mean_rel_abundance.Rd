% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adders_taxa.R
\name{add_mean_rel_abundance}
\alias{add_mean_rel_abundance}
\title{Add average relative abundances to taxa table}
\usage{
add_mean_rel_abundance(ta, condition = NULL, test = NULL)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{condition}{A condition variable (character).}

\item{test}{Differential abundance test to perform.}
}
\value{
A tidytacos object
}
\description{
\code{add_mean_rel_abundance()} adds mean relative abundance values
for each taxon to the taxa table, overall or per sample group.
}
\details{
If \code{condition} is specified, the mean relative abundances will be calculated
separately for each group defined by the condition variable. This variable
should be present in the sample table.

If \code{condition} is specified, differential abundance testing can be performed
by setting the \code{test} argument. Options are \code{NULL} (default), \code{"wilcox"} or
\code{"t-test"}.
}
\seealso{
Other taxa-modifiers: 
\code{\link{add_prevalence}()},
\code{\link{add_taxon_name}()},
\code{\link{add_taxon_name_color}()}
}
\concept{taxa-modifiers}
