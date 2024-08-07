x <- matrix(
  c(1500, 1300, 280, 356),
  ncol = 2
)
rownames(x) <- c("taxon1", "taxon2")
colnames(x) <- c("sample1", "sample2")
test_data <- create_tidytacos(x, 
    taxa_are_columns=FALSE)

test_that("Can add sample tibble to ta object",{
    sample <- c("sample1","sample2")
    environment <- c("food fermentation", "human stool")
    smp_tibble <- tibble::tibble(sample, environment)
    suppressMessages(
        test_data <- test_data %>% add_metadata(smp_tibble)
    )
    expect_true("environment" %in% names(test_data$samples))
})

test_that("Can add lib sizes", {
    ta_lib <- test_data %>% add_total_count()
    expect_equal(ta_lib$samples$total_count, c(2800, 636))
    
})

test_that("Can add alpha diversity metrics", {
    ta_alpha <- test_data %>% add_alphas()
    expect_equal(ta_alpha$samples$invsimpson, 
                    c(1.989, 1.971), tolerance=1e-3, ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$shannon, 
                    c(0.690, 0.686), tolerance=1e-3, ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$simpson, 
                    c(0.497, 0.493), tolerance=1e-3, ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$obs, 
                    c(2, 2), ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$s.chao1, 
                    c(2,2), ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$s.ace, 
                    c(NaN,NaN), ignore_attr=TRUE)
    expect_equal(ta_alpha$samples$pielou, 
                    c(0.996,0.990), tolerance=1e-3, ignore_attr=TRUE)
})

test_that("Can add spike-ratio", {

    ta_spike_ratio <- test_data %>% add_spike_ratio("t1")
    expect_true("spike_ratio" %in% names(ta_spike_ratio$samples))
})

test_that("Can perform anosim test", {

    suppressWarnings(
        anosim <- urt %>% perform_anosim("plate")
        ) # empty samples
    expect_gt(anosim$signif, 0.05)
})
