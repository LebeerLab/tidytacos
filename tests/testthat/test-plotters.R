# Tests could be expanded with proto testing of diff layers,
# but for now a heuristic test
# to see if the same plot is generate will suffice

test_that("Barplot returns identical plot", {
    bp  <- urt %>% tacoplot_stack()
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Default barplot", bp)
})

test_that("Barplot raises warning when aggregating samples", {
    expect_warning(bp <- urt %>% tacoplot_stack(n=5, x=participant))
    skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger("Custom barplot", bp)
})

test_that("Barplot raises error when providing non-existant label", {
    expect_error(urt %>% tacoplot_stack(x=imagined))
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
        urt %>% tacoplot_ord_ly(x=location)
    )

})

test_that("Can create venndiagram", {
    skip_if_not_installed("ggVenDiagram")
    venn <- urt %>% tacoplot_venn(location)
    vdiffr::expect_doppelganger("Venndiagram", venn)
})
