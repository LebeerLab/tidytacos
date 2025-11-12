# Return an euler diagram of overlapping taxon_ids between conditions

Return an euler diagram of overlapping taxon_ids between conditions

## Usage

``` r
tacoplot_euler(ta, condition, shape = "ellipse", ...)
```

## Arguments

- ta:

  A tidytacos object.

- condition:

  The name of a variable in the samples table that contains a
  categorical value.

- shape:

  shape to plot the groups in; choice from circle or ellipse

- ...:

  Arguments passed on to
  [`eulerr::euler`](https://jolars.github.io/eulerr/reference/euler.html)

  `combinations`

  :   set relationships as a named numeric vector, matrix, or data.frame
      (see **methods (by class)**)
