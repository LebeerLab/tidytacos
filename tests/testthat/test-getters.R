len_taxa <- length(urt$taxa$taxon_id)
len_samples <- 214 #some empty samples in the df
cols_samples <- 9

test_that("Correct numbers reported", {
    expect_snapshot(urt %>% tacosum())
})

test_that("Relative abundance matrix can be created", {
    relabun_mat <- urt %>% rel_abundance_matrix()
    expect_lte(max(relabun_mat), 1)
    expect_equal(dim(relabun_mat), c(len_samples, len_taxa))
})

test_that("Beta diversity table can be generated", {
    beta_table <- urt %>% betas()
    expect_lte(max(beta_table$beta), 1)
    expect_equal(dim(beta_table), c(22791, 1+cols_samples*2))
})

test_that("Occurences are calculated",{
    occ <- urt %>% occurrences("location")
    expect_equal(dim(occ), c(2*len_taxa, 3))
    expect_true(all(c("location","occurrence") %in% names(occ)))
    expect_equal(sum(occ$occurrence), 7693)
})

test_that("Presence/absense is calculated", {
    pres_ab <- urt %>% occurrences("location", pres_abs=TRUE)
    expect_equal(dim(pres_ab), c(4*len_taxa, 4))
    expect_true(all(c("location", "presence", "n") %in% names(pres_ab)))
    expect_equal(sum(pres_ab$n), len_samples * len_taxa)
})

test_that("Mean relative abundance is computed", {
    mra <- urt %>% mean_rel_abundances()
    expect_equal(sum(mra$mean_rel_abundance), 1)
    expect_equal(length(mra$taxon_id), len_taxa)
})

test_that("Mean relative abundance is computed for a condition", {
    mra <- urt %>% mean_rel_abundances("location")
    cond_opt <- unique(urt$samples$location)
    expect_equal(sum(mra$mean_rel_abundance), length(cond_opt))
    expect_equal(length(mra$taxon_id), len_taxa*length(cond_opt))
})

test_that("Extract sample tibble", {
    expect_equal(urt %>% samples(), urt$samples)
})

test_that("Extract taxon tibble", {
    expect_equal(urt %>% taxa(), urt$taxa)
})

test_that("Extract abundances tibble", {
    expect_equal(urt %>% counts(), urt$counts)
})

test_that("Extract all three tables in one large tibble",{
    wide_ta <- urt %>% everything()
    expect_equal(dim(wide_ta), c(7693, 10+cols_samples))
})

test_that("Perform adonis shows stable output", {
    set.seed(42)
    adn <- urt %>% perform_adonis(c("location", "method"))
    expect_snapshot(adn)
})

test_that("Can create a list of unique taxa per condition", {
    taxa_list <- urt %>% taxonlist_per_condition(location)
    expect_snapshot(taxa_list)
})
