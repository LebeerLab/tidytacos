% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/converters.R
\name{create_tidytacos}
\alias{create_tidytacos}
\title{Initiate tidytacos object}
\usage{
create_tidytacos(counts_matrix, taxa_are_columns = TRUE)
}
\arguments{
\item{counts_matrix}{Numerical matrix containing the count data.}

\item{taxa_are_columns}{A logical scalar. Are the taxa defined in columns?}
}
\description{
\code{tidytacos} returns a tidytacos object given a numeric matrix.
}
\details{
This function initiates a tidytacos object based on a numeric matrix. It
will automatically create a dummy taxa table and sample table which will need to be
updated using the function \code{\link{add_metadata}}.
}
\examples{
# Initiate count matrix
x <- matrix(
 c(1500, 1300, 280, 356),
 ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x,
                     taxa_are_columns = FALSE)


\dontrun{
tidytacos("a")
}

}
