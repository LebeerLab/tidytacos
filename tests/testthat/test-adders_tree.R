
test_that("Can add tree to tidytacos object", {

    urt_tree <- urt %>% 
      filter_samples(sample=="CON9-NF-S") %>% 
      add_tree()   
    
    expect_true(is(urt_tree$tree, "phylo"))
})
