% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lda_models.R
\name{perform_lda}
\alias{perform_lda}
\title{LDA model estimation}
\usage{
perform_lda(
  ta,
  k,
  min_prevalence = 0.05,
  taxon = taxon_id,
  sample = sample_id,
  ...
)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{k}{The number of topics to estimate.}

\item{min_prevalence}{The lowest percentage (0-1) of samples
taxa need to be present in for the taxa to be used in model estimation.}

\item{taxon}{The column name in the taxa table with taxa identifiers.}

\item{sample}{The column name in the sample table with sample identifiers.}

\item{...}{
  Arguments passed on to \code{\link[topicmodels:lda]{topicmodels::LDA}}
  \describe{
    \item{\code{x}}{Object of class \code{"DocumentTermMatrix"} with
    term-frequency weighting or an object coercible to a
    \code{"simple_triplet_matrix"} with integer entries.}
    \item{\code{method}}{The method to be used for fitting; currently 
      \code{method = "VEM"} or \code{method= "Gibbs"} are
      supported.}
    \item{\code{control}}{A named list of the control parameters for estimation
    or an object of class \code{"LDAcontrol"}.}
    \item{\code{model}}{Object of class \code{"LDA"} for initialization.}
  }}
}
\value{
A list of estimated topics, terms and the perplexity of the model.
}
\description{
\code{perform_lda()} estimates a Latent Dirichlet Allocation (LDA) model
on the counts matrix of a tidytacos object. The function returns
the estimated topics and terms, as well as the perplexity of the model.
}
