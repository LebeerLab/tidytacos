# RAREFY
test_that("Raise error on rarification higher
than lowest abundance when replace=FALSE.", {
expect_error(urt %>% rarefy(100))
})

test_that("Max abundance equals n", {
    ab <- urt %>% rarefy(10) %>%
    counts() %>% dplyr::pull(count)
    expect_equal(max(ab), 10)
})

test_that("Max abundance equals n with replace=TRUE", {
    ab <- urt %>% rarefy(100, replace=TRUE) %>%
    counts() %>% dplyr::pull(count)
    expect_equal(max(ab), 100)
})

# CHANGE_ID_SAMPLES
test_that("Change sample ID based on unique column",{
    uq_part <- urt %>% filter_samples(location == "NF", method=="A")
    # remove empty samples, as they will cause a mismatch
    uq_part <- uq_part %>% remove_empty_samples()
    uq_part_id <- uq_part %>% change_id_samples(participant)
    expect_true(identical(
        sort(unique(uq_part$samples$participant)),
        sort(unique(uq_part_id$counts$sample_id))
    ))
})

test_that("Raise error on non-unique sample_id",{
    expect_error(urt %>% change_id_samples(condition))
})

# CHANGE_ID_TAXA

test_that("Change sample ID based on unique column",{

    urt_id <- urt %>% change_id_taxa(sequence)
    expect_true(identical(
        sort(unique(urt_id$counts$taxon_id)),
        sort(unique(urt$taxa$sequence))
    ))
})

test_that("Raise error on non-unique taxon_id",{
    expect_error(urt %>% change_id_taxa(genus))
})

# AGG SAMPLES
nsamples <- 217
nsamples_agg <- 139

test_that("No aggregation of sample table when
    metadata is unique for every row", {
    sample_id_before <- urt$samples$sample_id
    ta_agg <- urt %>% aggregate_samples()
    sample_id_after <- ta_agg$samples$sample_id
    expect_length(sample_id_after, nsamples)
    expect_length(sample_id_before, nsamples)
})

test_that("Aggregation of sample table succeeds", {
    sample_id_before <- urt$samples$sample_id
    # Make one column not unique to test agg
    ta_agg <- urt %>%
        mutate_samples(location = "NF") %>%
        # sample name is unique, so remove for this test
        select_samples(-sample) %>%
        aggregate_samples()
    sample_id_after <- ta_agg$samples$sample_id
    expect_length(sample_id_before, nsamples)
    expect_length(sample_id_after, nsamples_agg)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta_agg$counts$sample_id)),
        nsamples_agg
    )
})

# AGG TAXA
ntaxa <- 1957
ntaxa_agg <- 1066
ntaxa_agg_phylum <- 32

test_that("No aggregation of taxa table when metadata
    is unique for every row", {
    taxon_id_before <- urt$taxa$taxon_id
    ta_agg <- urt %>% aggregate_taxa()
    taxon_id_after <- ta_agg$taxa$taxon_id
    expect_length(taxon_id_after, ntaxa)
    expect_length(taxon_id_before, ntaxa)
})

test_that("Aggregation of taxa table succeeds when no rank is specified", {
    taxon_id_before <- urt$taxa$taxon_id
    # Make one column not unique to test agg
    ta_agg <- urt %>%
        mutate_taxa(sequence = 0) %>%
        aggregate_taxa()
    taxon_id_after <- ta_agg$taxa$taxon_id
    expect_length(taxon_id_before, ntaxa)
    expect_length(taxon_id_after, ntaxa_agg)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta_agg$counts$taxon_id)),
        ntaxa_agg
    )
})

test_that("Aggregation of taxa table succeeds when rank is specified", {
    taxon_id_before <- urt$taxa$taxon_id
    # Make one column not unique to test agg
    ta_agg <- urt %>% aggregate_taxa(rank="phylum")
    taxon_id_after <- ta_agg$taxa$taxon_id
    expect_length(taxon_id_before, ntaxa)
    expect_length(taxon_id_after, ntaxa_agg_phylum)
    # test abundances table for irregularities
    expect_lte(
        length(unique(ta_agg$counts$taxon_id)),
        ntaxa_agg_phylum
    )
})

# TRIM ASVS

test_that("Trim ASVs works", {
    ta_trim <- urt %>% trim_asvs(1, 5)
    longest_seq <- max(str_length(ta_trim$taxa$sequence))
    expect_lte(longest_seq, 5)
})


# SELECT SAMPLES

test_that("Removing sample_id throws error", {
expect_error(urt %>% select_samples(-sample_id))
})

test_that("Removing a column from samples is succesful", {
    expect_true("method" %in% names(urt$samples))
    urt_t <- urt %>% select_samples(-method)
    expect_true(!"method" %in% names(urt_t$samples))
})

# SELECT TAXA

test_that("Removing taxon_id throws error", {
expect_error(urt %>% select_taxa(-taxon_id))
})

test_that("Removing a column from taxa is succesful", {
    expect_true("phylum" %in% names(urt$taxa))
    urt_t <- urt %>% select_taxa(-phylum)
    expect_true(!"phylum" %in% names(urt_t$taxa))
})

# SELECT ABUNDANCES

# Add removable column
urt$counts <- urt$counts %>% mutate(remove_me = 0)

test_that("Removing sample_id throws error", {
expect_error(urt %>% select_counts(-sample_id))
})

test_that("Removing taxon_id throws error", {
expect_error(urt %>% select_counts(-taxon_id))
})

test_that("Removing abundance column throws error", {
expect_error(urt %>% select_counts(-abundance))
})

test_that("Removing a column from abundances is succesful", {
    expect_true("remove_me" %in% names(urt$counts))
    urt_t <- urt %>% select_counts(-remove_me)
    expect_true(!"remove_me" %in% names(urt_t$counts))
})

# MUTATE TAXA
test_that("Mutating taxa", {
    ta_gs <- urt %>% mutate_taxa(new_col = paste(genus, species))
    expect_true("new_col" %in% names(ta_gs$taxa))
})

test_that("Mutating taxa, where taxon_id is removed throws error", {
    expect_error(urt %>% mutate_taxa(taxon_id = NULL))
})

# MUTATE SAMPLES
test_that("Mutating samples", {
    ta_con <- urt %>% mutate_samples(new_col = paste(condition, location))
    expect_true("new_col" %in% names(ta_con$samples))
})

test_that("Mutating samples, where sample_id is removed throws error", {
    expect_error(urt %>% mutate_samples(sample_id = NULL))
})

# MUTATE SAMPLES
test_that("Mutating abundances", {
    ta_st <- urt %>% mutate_counts(new_col = paste(sample_id, taxon_id))
    expect_true("new_col" %in% names(ta_st$counts))
})

test_that("Mutating abundances, where sample_id is removed throws error", {
    expect_error(urt %>% mutate_counts(sample_id = NULL))
})

test_that("Mutating abundances, where taxon_id is removed throws error", {
    expect_error(urt %>% mutate_counts(taxon_id = NULL))
})

test_that("Mutating abundances, where abundances is removed throws error", {
    expect_error(urt %>% mutate_counts(count = NULL))
})

test_that("Filtering out all samples raises error", {
    expect_error(urt %>% filter_samples(condition == "Non-existent"))
})

test_that("Filtering returns the expected amount of samples", {
    ta_f <- urt %>% filter_samples(method == "S", location == "NF")
    expect_equal(length(ta_f$samples$sample_id), 96)
})

test_that("Filtering out all taxa raises error", {
    expect_error(urt %>% filter_taxa(species == "Impossiblum"))
})

test_that("Filtering returns the expected amount of taxa", {
    ta_f <- urt %>% filter_taxa(species == "coli")
    expect_equal(length(ta_f$taxa$taxon_id), 1)
})

test_that("Filtering out all abundances raises error", {
    expect_error(urt %>% filter_counts(count == "Impossible"))
})

test_that("Filtering returns the expected amount of samples", {
    ta_f <- urt %>% filter_counts(count >= 10000)
    expect_equal(length(ta_f$counts$sample_id), 70)
    expect_equal(length(ta_f$taxa$taxon_id), 18)
})
