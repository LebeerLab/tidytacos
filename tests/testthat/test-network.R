# SparCC - MCL
test_that("Can run sparcc-MCL worfklow", {
skip_if_not_installed("MCL")
skip_if_not_installed("SpiecEasi")
skip_if_not_installed("Matrix")
urt_clust <- urt %>% add_sparcc_network_clusters(taxon_name=sequence)
expect_true("cluster" %in% names(urt_clust$taxa))

})

# Eigentaxa
test_that("Can add eigentaxa", {
    skip_if_not_installed("MCL")
    skip_if_not_installed("SpiecEasi")
    skip_if_not_installed("Matrix")
    urt_eigen <- urt %>% add_eigentaxa(taxon_name=sequence)
    expect_gt(length(urt_eigen$samples), length(urt$samples))
})
