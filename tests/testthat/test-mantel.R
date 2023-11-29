test_that("Can perform Mantel test on one variable", {

    expect_warning(
        mantel <- urt %>% perform_mantel_test("plate")
    ) # empty samples
    expect_gt(mantel$signif, 0.2)

})

test_that("Can perform Mantel test on array of variables", {

    expect_warning(
        mantel <- urt %>% perform_mantel_test(c("location","method"))
    ) # empty samples
    expect_gt(mantel$signif, 0.2)

})


test_that("Throw error when trying to test one constant variable", {

    expect_error(urt %>% perform_mantel_test("run"))
    
})
