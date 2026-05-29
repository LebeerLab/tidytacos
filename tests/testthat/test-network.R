# network
test_that("Can run network method with p-val calculation", {

  uns <- urt %>% 
  filter_samples(location=="NF", method=="S")

  uns %>%
  network(taxon_name=taxon_name, calculate_p=TRUE, Nboot=2) %>%
  expect_warning("taxon_name not found in the taxa table, resorting to add_taxon_name") %>%
  expect_warning("Some p-values are NA, consider increasing the number of bootstraps")  

})

# SparCC - MCL
test_that("Can run sparcc-MCL worfklow", {
skip_if_not_installed("MCL")
skip_if_not_installed("SpiecEasi")
skip_if_not_installed("Matrix")
urt_clust <- cluster_taxa(urt, taxon_name = sequence)
expect_true("cluster" %in% names(urt_clust$taxa))

})

# Eigentaxa
test_that("Can add eigentaxa", {
  skip_if_not_installed("MCL")
  skip_if_not_installed("SpiecEasi")
  skip_if_not_installed("Matrix")
  urt_eigen <- add_eigentaxa(urt, taxon_name = sequence)
  expect_gt(length(urt_eigen$samples), length(urt$samples))
})
