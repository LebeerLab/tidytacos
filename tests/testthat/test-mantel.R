test_that("Can perform Mantel test on one variable", {

    suppressWarnings(
        mantel <- urt %>% perform_mantel_test("plate")
    ) # empty samples
    expect_gt(mantel$signif, 0.1)

})

test_that("Can perform Mantel test on array of variables", {

    skip_if_not_installed("fastDummies")
    expect_warning(
        expect_warning(
            expect_warning(
                mantel <- urt %>% 
                  perform_mantel_test(c("location","method")),
                regexp = "Empty samples found"
            ),
        regexp = "Removed 3 empty samples"
        ),
        regexp = "non-square matrix"
    )
    expect_gt(mantel$signif, 0.1)

})


test_that("Throw error when trying to test one constant variable", {

    expect_error(
        suppressWarnings(
            urt %>% perform_mantel_test("run")
        )
    )
    
})
