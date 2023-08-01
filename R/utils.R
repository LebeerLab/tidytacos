#' Create a tidytacos object for testing/example purporses
#' @export
#' @return a small tidytacos object
create_test_ta <- function(){
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
#' @param ta a tidytacos object
#' @export
#' @return the tidytacos object minus the empty samples
remove_empty_samples <- function(ta){
  present_samples <- ta %>% abundances() %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::count() %>% pull(sample_id)
  empty_samples <- ta$samples$sample_id[!(ta$samples$sample_id %in% present_samples)] 
  ta <- ta %>% filter_samples(!sample_id %in% empty_samples)
  ta
}