# Tests could be expanded with proto testing of diff layers,
# but for now a heuristic test
# to see if the same plot is generate will suffice

## TACOPLOT_STACK
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
        suppressWarnings(length(urt %>% remove_empty_samples() %>% samples %>% pull(sample)))  
    )
})

test_that("Barplot raises warning when aggregating samples", {
    testthat::skip_on_ci()
    expect_warning(bp <- urt %>% tacoplot_stack(n=5, x=participant))
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Custom barplot", bp)
})

test_that("Barplot raises error when providing non-existant label", {
    expect_error(
        urt %>%
        tacoplot_stack(x = imagined)
    )
})

test_that("Pieplot can visualize 1 sample", {
    
    pp <- urt %>%
    filter_samples(sample_id == "s1") %>%
    tacoplot_stack(pie = TRUE)
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Pieplot", pp)
})

test_that("Pieplot raises error when visualizing more than 1 sample", {
    expect_error(urt %>% tacoplot_stack(pie=TRUE))
})

test_that("Can run tacoplot zoom on a sample",{
    tzoom <- urt %>% 
      filter_samples(sample_id == "s3") %>% 
      tacoplot_zoom()
    vdiffr::expect_doppelganger("Zoomed barplot", tzoom)
})


## TACOPLOT_ORD
test_that("get_ord_stat returns correct stats",{
    stat <- urt %>% get_ord_stat("location", stat.method="permanova")
    expect_equal(stat$signif, 0.001)
})

test_that("Tacoplot_ord works with tsne", {
    skip_if_not_installed("Rtsne")
    expect_no_error(
        expect_warning(
            urt %>% tacoplot_ord(x=location, ord="tsne"),
            "Removed 3 empty samples."
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
            urt %>% tacoplot_ord_ly(x=location),
            "Removed 3 empty samples."
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

## TACOPLOT_VENN
test_that("Can create venndiagram", {
    skip_if_not_installed("ggVennDiagram")
    testthat::skip_on_ci()
    venn <- urt %>% tacoplot_venn(location)
    vdiffr::expect_doppelganger("Venndiagram", venn)
})

test_that("Can create interactive venndiagram",{
    skip_if_not_installed("ggVennDiagram")
    skip_if_not_installed("plotly")
    testthat::skip_on_ci()
    suppressWarnings(venn_ly <- urt %>% tacoplot_venn_ly(plate))
    vdiffr::expect_doppelganger("Interactive Venndiagram", venn_ly)
})

test_that("Can create euler plot", {
    skip_if_not_installed("eulerr")
    testthat::skip_on_ci()
    euler <- urt %>% tacoplot_euler(location)
    expect_equal(euler$name, "euler.diagram")
    expect_equal(rownames(euler$data$ellipses), c("NF","N"))
})

## TACOPLOT_ALPHAS
test_that("Tacoplot_alphas works", {
    expect_no_error(
        expect_warning(
            urt %>% tacoplot_alphas(location),
            "Removed 3 empty samples."
        )
    )
})

## TACOPLOT_PREVALENCES

test_that("Tacoplot_prevalences works", {
    skip_if_not_installed("rstatix")
    p <- urt %>% tacoplot_prevalences(location)
    expect_equal(class(p), "pheatmap")
})
