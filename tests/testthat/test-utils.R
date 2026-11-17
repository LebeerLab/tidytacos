x <- matrix(
   c(1500, 1300, 14, 280, 356, 9, 1500, 1300, 14),
   ncol = 3
)
rownames(x) <- c("taxon1", "taxon2", "taxon3")
colnames(x) <- c("sample1", "sample2", "sample1")

# Convert to tidytacos object
dupl_data <- create_tidytacos(x,
taxa_are_columns = FALSE
)


test_that("Can create test taco", {
    tt <- test_taco()
    expect_snapshot(tt)
})

test_that("Complain about missing packages", {
    expect_snapshot(force_optional_dependency("dadjokeapi"), error=T)
})

test_that("Can remove duplicate samples", {
    
   dedup <- remove_duplicate_samples(dupl_data)
   expect_equal(length(dedup$samples$sample_id), 2)
   expect_equal(length(dupl_data$samples$sample_id), 3)

})
