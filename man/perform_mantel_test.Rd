% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mantel.R
\name{perform_mantel_test}
\alias{perform_mantel_test}
\title{Determine the correlation between the distance of the counts in a tidytacos object and a sample variable, multiple sample variables or another matrix.}
\usage{
perform_mantel_test(ta, comparison, ...)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{comparison}{A distance to compare against. This can be any of the following:
\itemize{
\item The name of the variable in the sample table to use for comparison
\item A list of names of variables in the sample table.
\item A distance matrix.
}}

\item{...}{Additional arguments to pass to the mantel function.}
}
\value{
The mantel test statistics
}
\description{
This function performs a mantel test using the dissimilarity matrix
of the tidytacos object supplied and a second distance matrix generated from user input.
}
