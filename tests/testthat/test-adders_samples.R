x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")
test_data <- create_tidytacos(x, 
    taxa_are_columns=FALSE)

test_that("Can add sample tibble to ta object",{
    sample <- c("sample1","sample2")
    environment <- c("food fermentation", "human stool")
    smp_tibble <- tibble::tibble(sample, environment)
    suppressMessages(
        test_data <- test_data %>% add_metadata(smp_tibble)
    )
    expect_true("environment" %in% names(test_data$samples))
})

test_that("Can add lib sizes", {
    ta_lib <- test_data %>% add_total_count()
    expect_equal(ta_lib$samples$total_count, c(2800, 636))
    
})

test_that("Can add alpha diversity metrics", {
    ta_alpha <- test_data %>% add_alphas()
    expect_equal(ta_alpha$samples$invsimpson, 
                    c(1.989, 1.971), tolerance=1e-3, ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$shannon, 
                    c(0.690, 0.686), tolerance=1e-3, ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$simpson, 
                    c(0.497, 0.493), tolerance=1e-3, ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$obs, 
                    c(2, 2), ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$s.chao1, 
                    c(2,2), ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$s.ace, 
                    c(NaN,NaN), ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$pielou, 
                    c(0.996,0.990), tolerance=1e-3, ignore_attr=TRUE)
})

test_that("Can add alpha diversities with subsampling", {
    ta_a_sub <- test_data %>%
      add_alphas(subsample=TRUE, min_lib_size=1000)
    expect_equal(ta_a_sub$samples$mean_shannon, .690, tolerance=.001)
    expect_equal(ta_a_sub$samples$median_obs, 2, tolerance=.001)
})

test_that("Can add alpha diversity of a single sample", {
 
  single_sample <- test_data %>%
    filter_samples(sample_id == "s1")
  inv_simp <- single_sample %>%
    add_alpha() %>% 
    samples() %>% 
    pull("invsimpson")
  expect_equal(inv_simp, 1.99, tolerance = .01)
})

test_that("Can add spike-ratio", {

    ta_spike_ratio <- test_data %>% add_spike_ratio("t1")
    expect_true("spike_ratio" %in% names(ta_spike_ratio$samples))
})

test_that("Can add dominant taxa to the table", {
    urt_dom <- urt %>% add_dominant_taxa()
    expect_equal(sum(is.na(urt_dom$samples$dominant_taxon)), 149)
    # try again with other threshold using previous ttaco
    # this test ensures that no new variables like dominant_taxon.x start popping up

    urt_dom2 <- urt_dom %>% add_dominant_taxa(threshold_dominance=.3, taxon_name="phylum")
    expect_equal(sum(is.na(urt_dom2$samples$dominant_taxon)), 73)

})


test_that("Add dominant taxa errors when trying to use nonexistent taxon column", {
    expect_error(
      urt_dom <- urt %>% add_dominant_taxa(taxon_name="idunno"),
      "Taxa table requires a column idunno that defines the taxon name."
    )  
})

test_that("Can perform anosim test", {

    suppressWarnings(
        anosim <- urt %>% perform_anosim("plate")
        ) # empty samples
    expect_gt(anosim$signif, 0.05)
})
