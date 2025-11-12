# Construct a phylogeny from ASV sequences

`add_tree()` performs a multiple sequence alignment on the ASV sequences
in the taxa table and then constructs a phylogenetic tree. The tree is
added to the tidytacos object under a variable called `tree`.

## Usage

``` r
add_tree(
  ta,
  sequence_var = sequence,
  aln_args = list(),
  tree_fit_args = list()
)
```

## Arguments

- ta:

  A tidytacos object.

- sequence_var:

  The name of the column in the taxa table that contains the ASV
  sequences

- aln_args:

  An optional list of arguments to pass to the `DECIPHER::AlignSeqs()`
  function

- tree_fit_args:

  An optional list of arguments to pass to the
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  function

## Value

A tidytacos object with a phylotree slot

## See also

Other unifrac-distance-functions:
[`calculate_unifrac_distances()`](https://lebeerlab.github.io/tidytacos/reference/calculate_unifrac_distances.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# filter taxa to speed up calculation time
urt_sub <- urt %>% filter_taxa(taxon_id %in% head(urt$taxa$taxon_id))
# infer tree
urt_sub_tree <- add_tree(urt_sub,
                         sequence_var="sequence",
                         aln_args=list(verbose=FALSE)
                        )
# inspect tree
urt_sub_tree$tree
} # }
```
