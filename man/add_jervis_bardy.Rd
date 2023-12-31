% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adders_taxa.R
\name{add_jervis_bardy}
\alias{add_jervis_bardy}
\title{Apply the taxon QC method of Jervis-Bardy}
\usage{
add_jervis_bardy(ta, dna_conc, sample_condition = T, min_pres = 3)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{dna_conc}{A variable in the samples table that contains dna
concetrations (unquoted).}

\item{sample_condition}{An optional extra condition that samples must pass
before calculations.}

\item{min_pres}{The minimum number of samples a taxon has to be present in
for its correlation to be calculated.}
}
\description{
\code{add_jervis_bardy} calculates the spearman correlation between relative abundance and
sample DNA concentration, for each taxon and adds the correlation metric and p-value to the taxa table under the column names "jb_cor" and "jb_p", respectively. If taxa show a distribution that is negatively correlated with DNA concentration, it indicates their potential as contaminants.
}
\details{
See:
J. Jervis-Bardy et al., “Deriving accurate microbiota profiles from
human samples with low bacterial content through post-sequencing processing
of Illumina MiSeq data,” Microbiome, vol. 3, no. 1, Art. no. 1, 2015, doi:
10.1186/s40168-015-0083-8.
}
