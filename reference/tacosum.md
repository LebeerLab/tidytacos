# Return some descriptive numbers

`tacosum()` returns the number of samples, taxa and reads in the
tidytacos object.

## Usage

``` r
tacosum(ta)
```

## Arguments

- ta:

  A tidytacos object.

## Examples

``` r
# Initiate counts matrix
x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x, taxa_are_columns = FALSE)

# Report numbers
numbers <- data %>% tacosum()
```
