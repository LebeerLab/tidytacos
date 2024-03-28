test_that("Can plot ord after taxa aggregation", {
    # Previously the issue was that the taxon field was removed during aggregation
    # and the plot_ord function uses this to get the relative abundance
    urta <- urt %>% aggregate_taxa(rank="family") 
    expect_no_error(
        urta %>% tacoplot_ord(x="plate")
    )
})