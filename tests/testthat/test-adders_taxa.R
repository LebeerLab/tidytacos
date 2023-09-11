x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")
test_data <- create_tidytacos(x,
    taxa_are_columns=FALSE)

test_that("Can reclassify using dada and small test database", {
    skip_if_not_installed("dada2")
    test_db = test_path("data/test_db.fa")
    expect_no_error(
        ta_reclass <- urt %>% classify_taxa(test_db)
    )
})

test_that("Can add taxon tibble", {
    taxon <- c("taxon1", "taxon2")
    genus <- c("Salmonella", "Lactobacillus")
    taxon_tibble <- tibble::tibble(taxon, genus)
    suppressMessages(
        test_data <- test_data %>% add_taxon_tibble(taxon_tibble)
    )
    expect_equal(test_data$taxa$genus, genus)
})

test_that("Can add taxon names", {
    expect_no_error(urt %>% add_taxon_name())
    expect_no_error(urt %>% add_taxon_name(include_species=TRUE))
    # TODO: find some more clever validation here
})

test_that("Add taxon name raises error when uring none existant method", {
    expect_error(urt %>% add_taxon_name(method="do something"))
    expect_error(urt %>% add_taxon_name(method="do something", include_species=TRUE))
})

test_that("Can add taxon name colors", {
    ta_col <- urt %>% add_taxon_name_color()
    expect_equal(length(levels(ta_col$taxa$taxon_name_color)),12)
})

test_that("Can add taxon name colors for specific number", {
    ta_col <- urt %>% add_taxon_name_color(n=4)
    expect_equal(length(levels(ta_col$taxa$taxon_name_color)),4)
})

test_that("Can add taxon name colors for specific taxa", {
    ta_col <- urt %>% add_taxon_name_color(taxa=c("t1","t2"))
    expect_equal(levels(ta_col$taxa$taxon_name_color),c("residual","Moraxella 1","Staphylococcus 1"))
})

test_that("Add taxon name color raises error when uring none existant method", {
    expect_error(urt %>% add_taxon_name_color(method="do something"))
    expect_error(urt %>% add_taxon_name_color(method="do something", include_species=TRUE))
})

test_that("Occurence can not be higher than amount of samples", {
    ta_occ <- urt %>% add_occurrences()
    expect_lte(max(ta_occ$taxa$occurrence), dim(urt$samples)[1])
})

test_that("Occurence in conditions with fischer test can be run", {
    ta_occ <- urt %>% add_occurrences(condition="location", fischer_test=TRUE)
    expect_true(all(c("occurrence_in_N","occurrence_in_NF","fisher_p") %in% names(ta_occ$taxa)))
    expect_lte(max(ta_occ$taxa$occurrence_in_NF), dim(urt$samples)[1])
    expect_lte(max(ta_occ$taxa$occurrence_in_N), dim(urt$samples)[1])
})

test_that("Relative occurences in conditions with fischer test can be run", {
    ta_occ <- urt %>% add_occurrences(condition="location", fischer_test=TRUE)
    expect_true(all(c("occurrence_in_N","occurrence_in_NF","fisher_p") %in% names(ta_occ$taxa)))
    expect_lte(max(ta_occ$taxa$occurrence_in_NF), dim(urt$samples)[1])
    expect_lte(max(ta_occ$taxa$occurrence_in_N), dim(urt$samples)[1])
})

test_that("Can add mean rel abundance", {
    ta_mr <- urt %>% add_mean_rel_abundances()
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance), 1)
})

test_that("Can add mean rel abundance with condition", {
    ta_mr <- urt %>% add_mean_rel_abundances(condition="location")
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance_in_N), 1)
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance_in_NF), 1)
})

test_that("Can add mean rel abundance with condition and wilcox tests", {
    ta_mr <- urt %>% add_mean_rel_abundances(condition="location", test="wilcox")
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance_in_N), 1)
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance_in_NF), 1)
})

test_that("Can add mean rel abundance with condition and t-test", {
    ta_mr <- urt %>% add_mean_rel_abundances(condition="location", test="t-test")
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance_in_N), 1)
    expect_equal(sum(ta_mr$taxa$mean_rel_abundance_in_NF), 1)
})

test_that("Raise error when non defined test is used in add mean rel abundance", {
    expect_error(urt %>% add_mean_rel_abundances(condition="location", test="bogus-test"))
})
