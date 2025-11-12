# Calculate unifrac distance matrix from a tidytacos object with a rooted tree

Calculate unifrac distance matrix from a tidytacos object with a rooted
tree

## Usage

``` r
calculate_unifrac_distances(ta, ...)
```

## Arguments

- ta:

  A tidytacos object with a rooted tree in the "tree" slot.

- ...:

  Arguments passed on to
  [`phyloseq::UniFrac`](https://rdrr.io/pkg/phyloseq/man/UniFrac-methods.html)

  `physeq`

  :   (Required).
      [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html),
      containing at minimum a phylogenetic tree
      ([`phylo-class`](https://rdrr.io/pkg/phyloseq/man/phylo-class.html))
      and contingency table
      ([`otu_table-class`](https://rdrr.io/pkg/phyloseq/man/otu_table-class.html)).
      See examples below for coercions that might be necessary.

  `weighted`

  :   (Optional). Logical. Should use weighted-UniFrac calculation?
      Weighted-UniFrac takes into account the relative abundance of
      species/taxa shared between samples, whereas unweighted-UniFrac
      only considers presence/absence. Default is `FALSE`, meaning the
      unweighted-UniFrac distance is calculated for all pairs of
      samples.

  `normalized`

  :   (Optional). Logical. Should the output be normalized such that
      values range from 0 to 1 independent of branch length values?
      Default is `TRUE`. Note that (unweighted) `UniFrac` is always
      normalized by total branch-length, and so this value is ignored
      when `weighted == FALSE`.

  `parallel`

  :   (Optional). Logical. Should execute calculation in parallel, using
      multiple CPU cores simultaneously? This can dramatically hasten
      the computation time for this function. However, it also requires
      that the user has registered a parallel “backend” prior to calling
      this function. Default is `FALSE`. If FALSE, UniFrac will register
      a serial backend so that `foreach::%dopar%` does not throw a
      warning.

  `fast`

  :   (Optional). Logical. DEPRECATED. Do you want to use the “Fast
      UniFrac” algorithm? Implemented natively in the
      `phyloseq-package`. `TRUE` is now the only supported option. There
      should be no difference in the output between the two algorithms.
      Moreover, the original UniFrac algorithm only outperforms this
      implementation of fast-UniFrac if the datasets are so small
      (approximated by the value of `ntaxa(physeq) * nsamples(physeq)`)
      that the difference in time is inconsequential (less than 1
      second). In practice it does not appear that this parameter should
      have ever been set to `FALSE`, and therefore the original UniFrac
      implementation perhaps never should have been supported here. For
      legacy code support the option is now deprecated here (the
      implementation was an internal function, anyway) and the `fast`
      option will remain for one release cycle before being removed
      completely in order to avoid causing unsupported-argument errors.

## Value

A distance matrix.

## See also

Other unifrac-distance-functions:
[`add_tree()`](https://lebeerlab.github.io/tidytacos/reference/add_tree.md)
