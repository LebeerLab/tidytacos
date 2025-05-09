% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/converters.R
\name{to_fasta}
\alias{to_fasta}
\title{Write the sequences of the taxa table to a fasta file}
\usage{
to_fasta(ta, filename = "asvs.fasta", taxon_col = taxon_id, seq_col = sequence)
}
\arguments{
\item{ta}{A tidytacos object.}

\item{filename}{The name of the resulting biom table file,
defaults to 'asvs.fasta'.}

\item{taxon_col}{The name of the column in the taxa table which
is to be used as id for the sequences (taxon_id by default).}

\item{seq_col}{The name of the sequence column in the taxa table
(sequence by default).}
}
\description{
Uses the taxon_col and sequence_col columns to write the sequences
into a fasta file per taxon.
}
\seealso{
Other export-methods: 
\code{\link{as_phyloseq}()},
\code{\link{to_biom}()},
\code{\link{write_tidytacos}()}
}
\concept{export-methods}
