test_that("Can perform Mantel test on one variable", {

    mantel <- urt %>% perform_mantel_test("plate")
    expect_gt(mantel$signif, 0.2)

})

test_that("Can perform Mantel test on array of variables", {

    mantel <- urt %>% perform_mantel_test(c("location","method"))
    expect_gt(mantel$signif, 0.2)

})


test_that("Throw error when trying to test one constant variable", {

    expect_error(urt %>% perform_mantel_test("run"))
    
})
