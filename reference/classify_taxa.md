# (Re)classify amplicon sequences

This function requires the DADA2 package to be installed.

## Usage

``` r
classify_taxa(
  ta,
  refdb,
  taxa = rep(T, times = length(taxon_id)),
  ranks = "default",
  sequence_var = "sequence",
  multithread = T,
  min_boot = 50,
  n_ranks = 7
)
```

## Arguments

- ta:

  A tidytacos object.

- refdb:

  The path to a DADA2-compatible reference database.

- taxa:

  An expression that specifies which taxa to (re)classify.

- ranks:

  A vector that specifies which ranks to (re)classify.

- sequence_var:

  The (quoted) name of a variable within the taxa table that contains
  (representative) sequences of the taxa.

- multithread:

  A boolean indicating whether to use multiple threads.

- min_boot:

  The minimum bootstrap value for taxonomy assignment.

- n_ranks:

  The number of ranks present in the reference database.

## Value

An updated tidytacos object.

## Details

`classify_taxa()` will (re)classify either all or a subset of the taxa,
given that a variable is present in the taxon table that contains
(representative) sequences of the taxa.

Ranks can be supplied as a named integer vector, where the names
represent taxonomic ranks and the integers represent positions of these
ranks in the taxonomy strings present in the reference database. Ranks
can also be supplied as just a character vector with the rank names; in
that case, it is assumed that the database taxonomy string follows the
default order (domain, phylum, class, order, family, genus, species). If
no ranks are supplied, taxa will be (re)classified at all default ranks.

## Examples

``` r
# we create a mock database
x <- c(
">Level1;Level2;Level3;Level4;Level5;Level6;",
"ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCATAGGATAC", 
">Level1;Level2;Level3;Level4;Level5;",
"CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC"
)
# and write it to a file
write(x, file="tmp-db.fna")

urt_reclass <- urt %>%
# filter out samples to save time
filter_samples(sample_id %in% c("s1","s2")) %>%
classify_taxa(
  "tmp-db.fna", n_ranks = 6, # the mock database is used here
  ranks=c("kingdom","phylum", "class", "order", "family", "genus")
)
# remove the temp file
unlink("tmp-db.fna")
```
