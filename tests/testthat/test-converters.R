test_that("Can read tidytacos object.", {
  ta <- read_tidytacos(test_path("data/urt"))
  expect_equal(attr(ta, "class"), "tidytacos")
})

test_that("Can save a tidytacos object.", {
  expect_no_warning(write_tidytacos(urt, "test"))
  on.exit(unlink("test", recursive=TRUE), add=TRUE, after=FALSE)
})

test_that("Can convert to phyloseq object.", {
  skip_if_not_installed("phyloseq")
  ta_phylo <- as_phyloseq(urt, sample=sample_id, taxon=taxon_id)
  expect_true(class(ta_phylo) == "phyloseq")
})

test_that("Can convert abundances to abundances matrix", {
  ab_mat <- urt %>% counts() %>% as_counts_matrix()
  expected_width <- 1957
  expected_height <- 214
  width <- dim(ab_mat)[2]
  height <- dim(ab_mat)[1]
  expect_equal(width, expected_width)
  expect_equal(height, expected_height)
})

test_that("Can merge two tidytacos", {
  ta_merged <- merge_tidytacos(urt, urt)
  
  final_sample_id <- "s434"
  expect_equal(
    ta_merged$samples$sample_id[length(ta_merged$samples$sample_id)],
    final_sample_id
  )
})
