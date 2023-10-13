#' Create a tidytacos object for testing/example purporses
#'
#' @return a small tidytacos object
#'
#' @export
create_test_ta <- function() {

  # Initiate abundance matrix
  x <- matrix(
        c(1500, 1300, 280, 356, 456, 678),
        ncol = 3
        )
  rownames(x) <- c("taxon1", "taxon2")
  colnames(x) <- c("sample1", "sample2", "sample3")

  # Convert to tidytacos object
  data <- create_tidytacos(x,
            taxa_are_columns = FALSE
          )

  data

}

#' Removes empty samples from the tidytacos object
#'
#' @param ta a tidytacos object
#' @return the tidytacos object minus the empty samples
#'
#' @export
remove_empty_samples <- function(ta){
  present_samples <- ta %>% counts() %>%
    dplyr::group_by(sample_id) %>%
    dplyr::count() %>% pull(sample_id)
  empty_samples <- ta$samples$sample_id[!(ta$samples$sample_id %in% present_samples)]
  ta <- ta %>% filter_samples(!sample_id %in% empty_samples)
  ta
}

# Checks if optional dependency is loaded and stops code if not.
force_optional_dependency <- function(optional_pkg, instructions=NULL){
  if (!requireNamespace(optional_pkg, quietly = TRUE)) {
      stop(paste("The", optional_pkg, "package must be installed to use this function.", instructions))
  }
  NULL
}
