# Tests could be expanded with proto testing of diff layers,
# but for now a heuristic test
# to see if the same plot is generate will suffice

test_that("Barplot returns identical plot", {
    testthat::skip_on_ci()
    bp  <- urt %>% tacoplot_stack()
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Default barplot", bp)
})

# unit test for bug where sample_name column would cause samples to aggregate on x axis.
test_that("Barplot does not aggregate when a 'sample_name' column exists in sample table",{
    urt_w_sample_name <- urt
    urt_w_sample_name$samples$sample_name <- urt$samples$sample

    plt_w_sample_name <- urt_w_sample_name %>% tacoplot_stack()

    expect_length(
        ggplot_build(plt_w_sample_name)$layout$panel_params[[1]]$x$get_labels(),
        length(urt %>% remove_empty_samples() %>% samples %>% pull(sample))  
    )
})

test_that("Barplot raises warning when aggregating samples", {
    testthat::skip_on_ci()
    expect_warning(bp <- urt %>% tacoplot_stack(n=5, x=participant))
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Custom barplot", bp)
})

test_that("Barplot raises error when providing non-existant label", {
    expect_error(urt %>% tacoplot_stack(x=imagined))
})

test_that("Tacoplot_ord works with tsne", {
    skip_if_not_installed("Rtsne")
    expect_no_error(
        expect_warning(
            urt %>% tacoplot_ord(x=location, ord="tsne")
        ) # empty samples
    )
})

test_that("Tacoplot_stack_ly works", {
    skip_if_not_installed("plotly")
    expect_no_error(
        bply <- urt %>% tacoplot_stack_ly()
    )
})

test_that("Tacoplot_ord_ly works", {
    skip_if_not_installed("plotly")
    expect_no_error(
        expect_warning(
            urt %>% tacoplot_ord_ly(x=location)
        ) # empty samples
    )

})

test_that("Tacoplot_ord_ly works with umap and 3 dims", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("umap")
    expect_no_error(
        expect_warning(
            urt %>% tacoplot_ord_ly(x=location, ord="umap", dims=3)
        ) # empty samples
    )

})

test_that("Can create venndiagram", {
    skip_if_not_installed("ggVennDiagram")
    testthat::skip_on_ci()
    venn <- urt %>% tacoplot_venn(location)
    vdiffr::expect_doppelganger("Venndiagram", venn)
})

test_that("Tacoplot_alphas works", {
    expect_no_error(
            urt %>% tacoplot_alphas(location)
    )
})
