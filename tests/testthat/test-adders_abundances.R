
x <- matrix(
   c(1500, 1300, 14, 280, 356, 9),
   ncol = 2
)
rownames(x) <- c("taxon1", "taxon2", "taxon3")
colnames(x) <- c("sample1", "sample2")

# Convert to tidytacos object
data <- create_tidytacos(x,
taxa_are_columns = FALSE
)
data$samples$spike_added <- c(100, 150)
data$samples$grams_source <- c(4,5)


test_that("Can add absolute abundance",{
   absolutes_expected <- c(10661, 4602, 9240, 5851, 100, 148)
   tdata <- data %>%
   add_absolute_abundance(spike_taxon = "t3")
   
   expect_contains(tdata$counts$absolute_abundance, absolutes_expected)
})

test_that("Can add density",{
   densities_expected <- c(2665, 920, 2310, 1170, 25, 30)
   tdata <- data %>%
   add_density(spike_taxon = "t3", material_sampled = grams_source)
   
   expect_contains(tdata$counts$density, densities_expected)
})