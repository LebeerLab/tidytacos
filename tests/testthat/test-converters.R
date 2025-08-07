test_that("Can read tidytacos object.", {
  testthat::expect_warning(ta <- read_tidytacos(test_path("data/urt")))
  expect_equal(attr(ta, "class"), "tidytacos")
})

test_that("Can read tidyamplicons object and informs on conversion", {
  expect_warning(ta <- read_tidytacos(test_path("data/ta_urt")))
  expect_equal(attr(ta, "class"), "tidytacos")
})

test_that("Can save a tidytacos object.", {
  expect_no_warning(write_tidytacos(urt, "test"))
  on.exit(unlink("test", recursive=TRUE), add=TRUE, after=FALSE)
})

test_that("Complains if not a matrix", {
    expect_snapshot(urt$counts %>% create_tidytacos(), error=T)
})

test_that("Can convert to phyloseq object.", {
  skip_if_not_installed("phyloseq")
  ta_phylo <- as_phyloseq(urt, sample=sample_id, taxon=taxon_id)
  expect_true(class(ta_phylo) == "phyloseq")
})

test_that("Can convert abundances to abundances matrix", {
  ab_mat <- urt %>% counts() %>% counts_matrix()
  expected_width <- 1957
  expected_height <- 214
  width <- dim(ab_mat)[2]
  height <- dim(ab_mat)[1]
  expect_equal(width, expected_width)
  expect_equal(height, expected_height)
})

test_that("Can merge two tidytacos", {
  expect_warning(
  ta_merged <- merge_tidytacos(urt, urt)
  )
  final_sample_id <- "s434"
  expect_equal(
    ta_merged$samples$sample_id[length(ta_merged$samples$sample_id)],
    final_sample_id
  )
})

test_that("Can reset ids", {
  urt_plate1 <- urt %>% 
                filter_samples(plate == 1)
  expect_equal(urt_plate1$samples$sample_id[1], "s3")
  
  urt_reset <- urt_plate1 %>% 
               reset_ids()
  expect_equal(urt_reset$samples$sample_id[1], "s1")

})

test_that("Can convert counts matrix to tibble", {

  tib <- urt %>% counts_matrix() %>% counts_tidy()
  expect_true(dplyr::is.tbl(tib))
})

test_that("Can write sequences to fasta file", {
  skip_if_not_installed("seqinr")
  urt %>% 
   filter_taxa(taxon_id == "t1") %>% 
   to_fasta("test.fasta")
  f <- seqinr::read.fasta("test.fasta")
  expect_equal(seqinr::GC(f$t1), 0.51, tolerance=0.01)
  unlink("test.fasta")
})

test_that("Can convert to biom format",{
  expect_warning(urt %>% to_biom("test.biom"),
                  "Removed 3 empty samples.")
  expect_true(file.exists("test.biom"))
})

test_that("Converted biom format can be read correctly",{
  skip_if_not_installed("biomformat")
  urtb <- biomformat::read_biom("test.biom")
  expect_equal(urtb$shape, c(1957, 214))
})

test_that("Can convert to phyloseq object",{
  skip_if_not_installed("phyloseq")
  expect_warning(

  urte <- urt %>% 
    remove_empty_samples() %>% 
    reset_ids(),
    regexp = "Removed 3 empty samples."
  )
  urtph <- urte %>% as_phyloseq()

  expect_no_error(u <- urtph %>% from_phyloseq())
  
})

test_that("Can convert from dada2",{

  NREADS <- 2114
  NTAXA <- 4
  NSAMPLES <- 2
  seqtab <- readRDS(system.file("extdata", "dada2", "seqtab.rds", package='tidytacos'))
  taxa <- readRDS(system.file("extdata", "dada2", "taxa.rds", package='tidytacos'))
  ta <- from_dada(seqtab, taxa)
  
  summ <- ta %>% tacosum()
  expect_equal(as.vector(summ), c(NSAMPLES, NTAXA, NREADS))
})
