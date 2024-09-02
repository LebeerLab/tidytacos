test_that("Can plot ord after taxa aggregation", {
    # Previously the issue was that the taxon field was removed during aggregation
    # and the plot_ord function uses this to get the relative abundance
    urta <- urt %>% aggregate_taxa(rank="family") 
    expect_no_error(
        expect_warning(
            urta %>% tacoplot_ord(x="plate")
        )
    )
})

test_that('Tacoplot stack accepts strings or objects', {

    expect_no_error(
        urt %>% tacoplot_stack(x="sample")
    )
    expect_no_error(
        urt %>% tacoplot_stack(x=sample)
    )

})

test_that("Aggregation works when correcting rank names", {

    t <- test_taco()
    t$taxa$domain <- c("Dom","Dom")
    t$taxa$phylum <- c("Phyl","Phyl")
    t$taxa$class <- c("Class1","Class2")
    t <- t %>% set_rank_names(c("domain","phylum","class"))
    t <- t %>% aggregate_taxa(rank="class")
    expect_equal(nrow(t$taxa), 1)
})