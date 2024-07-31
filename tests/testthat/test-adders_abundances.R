
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
   
   testthat::expect_contains(tdata$counts$absolute_abundance, absolutes_expected)
})

test_that("Show warning for samples without spike taxa",{
   
   w <- testthat::capture_warnings(data %>%
   add_absolute_abundance(spike_taxon = "t4"))
   testthat::expect_match(w, "Sample without spike taxon detected: s1")
   testthat::expect_match(w, "Sample without spike taxon detected: s2")
})

test_that("Throw error if spike_taxon not found",{
   colname <- "fake_col"
   error_expected <- paste(
      "Sample table requires a column",
      colname,
      "that defines the quantity of spike added to the sample."
   )
   testthat::expect_error(
      data %>% add_absolute_abundance(spike_added = fake_col), 
      error_expected, fixed=TRUE)
})

test_that("Can add density",{
   densities_expected <- c(2665, 920, 2310, 1170, 25, 30)
   tdata <- data %>%
   add_density(spike_taxon = "t3", material_sampled = grams_source)
   
   testthat::expect_contains(tdata$counts$density, densities_expected)
})

test_that("Throw error if spike_added is missing",{
   
   colname <- "fake_col"
   expected_error <- paste(
      "Sample table requires a column",
      colname,
      "that defines the quantity of spike added to the sample."
   )

   testthat::expect_error(
   data %>%
   add_density(spike_taxon = "t3", 
               material_sampled = grams_source, 
               spike_added = fake_col),
   expected_error, fixed = TRUE
   )
})

test_that("Throw error if material_sampled is missing",{
   
   colname <- "fake_col"
   expected_error <- paste(
      "Sample table requires a column",
      colname,
      "that defines the quantity of sample used."
   )

   testthat::expect_error(
   data %>%
   add_density(spike_taxon = "t3", 
               material_sampled = fake_col),
   expected_error, fixed = TRUE
   )
})

test_that("Show warning if samples miss spike",{
   w <- testthat::capture_warnings(data %>%
   add_density(spike_taxon = "t4", material_sampled = grams_source))
   
   testthat::expect_contains(w, 
   "Sample without spike taxon detected: s1 \nSample without spike taxon detected: s2 \n")
})
