# Clusters samples into n clusters using partitioning around mediods (PAM)

`add_pam_clusters()` clusters the samples into n clusters and adds these
clusters to a new variable "cluster " in the sample table.

## Usage

``` r
add_pam_clusters(ta, n_clusters)
```

## Arguments

- ta:

  A tidytacos object.

- n_clusters:

  The number of desired clusters.

## Value

A tidytacos object with a pam column in the samples table.

## Details

This function calculates the weighted Unifrac distance between samples
followed by PAM clustering of samples. The user provides a number of
desired clusters which will be used to assign the samples to. A new
variable named "cluster " will be added to the samples tibble of a
tidytacos object defining to what cluster a sample belongs.
