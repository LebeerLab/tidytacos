% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adders_taxa.R
\name{add_taxon_name_color}
\alias{add_taxon_name_color}
\title{Add taxon color for visualization.}
\usage{
add_taxon_name_color(
  ta,
  method = "mean_rel_abundance",
  n = 12,
  samples = NULL,
  taxa = NULL,
  rank = NULL,
  threshold_dominance = NULL
)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{method}{The method by which to arrange the taxon names. Currently only
mean_rel_abundance or dominance.}

\item{n}{An integer denoting the amount of most abundant taxa to display.
Capacity at 12.}

\item{samples}{An optional vector of sample_id's of interest.}

\item{taxa}{An optional vector of taxon_id's of interest.}

\item{rank}{An optional rank to aggregate taxa on.}

\item{threshold_dominance}{An optional threshold for the dominance method.}
}
\value{
A tidytacos object.
}
\description{
\code{add_rel_abundance()} determines the most abundant taxa and assigns
them a color for consistent color codes of each taxon in visualizations.
A rank can be supplied to aggregate colors higher than the current rank.
}
\examples{
# display the 5 most abundant taxa at genus lvl
urt \%>\% add_taxon_name_color(n=5, rank='genus') \%>\% tacoplot_stack()
}
\seealso{
Other taxa-modifiers: 
\code{\link{add_mean_rel_abundance}()},
\code{\link{add_prevalence}()},
\code{\link{add_taxon_name}()}
}
\concept{taxa-modifiers}
