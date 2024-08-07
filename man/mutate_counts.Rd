% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handlers.R
\name{mutate_counts}
\alias{mutate_counts}
\title{Create extra variables in the count table}
\usage{
mutate_counts(ta, ...)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{...}{Mutate criteria for the counts table.}
}
\value{
A tidytacos object.
}
\description{
Create extra variables in the count table
}
\examples{
# add a column to the counts table
urt \%>\% mutate_counts(cl10 = log10(count))
}
