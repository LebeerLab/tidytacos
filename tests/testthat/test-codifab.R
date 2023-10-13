ta_codifab <- urt %>% add_codifab(condition=location)

test_that("Can add logratios.", {
  ta_log <- urt %>% add_logratios()
  expect_true(tibble::is_tibble(ta_log$logratios))
})

test_that("Can perform adonis via the condition argument", {
    expect_true(tibble::is_tibble(ta_codifab$taxon_pairs))
    expect_lte(max(ta_codifab$taxon_pairs$wilcox_p), 1)
    expect_false(is.infinite(max(ta_codifab$taxon_pairs$wilcox_p)))
    expect_true("NF_vs_N" %in% names(ta_codifab$taxon_pairs))
})

test_that("Can perform adonis via the conditions argument", {
    ta_codifab_c <- urt %>% 
        add_codifab(condition=plate, conditions = c("1","2"))
    expect_true(tibble::is_tibble(ta_codifab_c$taxon_pairs))
    expect_lte(max(ta_codifab_c$taxon_pairs$wilcox_p), 1)
    expect_false(is.infinite(max(ta_codifab_c$taxon_pairs$wilcox_p)))
    expect_true("1_vs_2" %in% names(ta_codifab_c$taxon_pairs))
})

test_that("Throws an error when the condition variable does not exist in the sample table", {
    expect_error(
        urt %>% add_codifab(condition=something_weird),
        "condition field does not exist in sample table",
        fixed=TRUE
    )
})

test_that("Throws error when one or both conditions can't be found", {
    expect_error(
        urt %>% add_codifab(condition=plate, conditions=c("that_one","this_one")),
        "one or both conditions not found",
        fixed=TRUE
    )
    expect_error(
        urt %>% add_codifab(condition=plate, conditions=c("1","this_one")),
        "one or both conditions not found",
        fixed=TRUE
    )
})

test_that("Throws error when there is multiple conditions in the selected condition", {
    expect_error(
        urt %>% add_codifab(condition=plate),
        "there need to be exactly two conditions",
        fixed=TRUE
    )
})

test_that("Throws error when plotting the codifab without taxon pair data", {
    expect_error(
        urt %>% codifab_plot(nope),
        "Please first run add_codifab() to generate the taxon pair comparisons.",
        fixed=TRUE
    )
})

test_that("Throws error when plotting the codifab without valid comparison field", {
    expect_error(
        ta_codifab %>% codifab_plot(nope),
        " is not an existing comparison in the taxon_pairs table.$"
    )
})

test_that("Can plot identical plot without any errors", {
    # vdiffr::expect_doppelganger("Codifab Plot", ta_codifab %>% codifab_plot(NF_vs_N))
    # not deterministic in order of taxa
    expect_no_error(ta_codifab %>% codifab_plot(NF_vs_N))
})
