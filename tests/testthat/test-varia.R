test_that("Gives default rank names when not provided", {
  expect_equal(urt %>% rank_names(), 
  c("kingdom", "phylum", "class", "order", "family", "genus")
  )
})

test_that("Gives custom rank names when provided", {
  urt_r <- urt %>% set_rank_names(
    c("domain","phylum","order","family", "genus", "species")
  ) 
  expect_equal(urt_r %>% rank_names(), 
  c("domain","phylum","order","family", "genus", "species")
  )
})