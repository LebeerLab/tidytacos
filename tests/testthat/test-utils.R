test_that("Can create test taco", {
    tt <- test_taco()
    expect_snapshot(tt)
})

test_that("Complain about missing packages", {
    expect_snapshot(force_optional_dependency("dadjokeapi"), error=T)
})