# Write the sequences of the taxa table to a fasta file

Uses the taxon_col and sequence_col columns to write the sequences into
a fasta file per taxon.

## Usage

``` r
to_fasta(ta, filename = "asvs.fasta", taxon_col = taxon_id, seq_col = sequence)
```

## Arguments

- ta:

  A tidytacos object.

- filename:

  The name of the resulting biom table file, defaults to 'asvs.fasta'.

- taxon_col:

  The name of the column in the taxa table which is to be used as id for
  the sequences (taxon_id by default).

- seq_col:

  The name of the sequence column in the taxa table (sequence by
  default).

## See also

Other export-methods:
[`as_phyloseq()`](https://lebeerlab.github.io/tidytacos/reference/as_phyloseq.md),
[`to_biom()`](https://lebeerlab.github.io/tidytacos/reference/to_biom.md),
[`write_tidytacos()`](https://lebeerlab.github.io/tidytacos/reference/write_tidytacos.md)
