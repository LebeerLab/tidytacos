# Apply a sample filtering to the taxon and counts tables
process_sample_selection <- function(ta) {

  # filter counts table
  selected_samples <- ta$samples$sample_id
  ta$counts <- ta$counts %>%
    filter(sample_id %in% selected_samples)

  # filter taxon table
  selected_taxa <- ta$counts$taxon_id %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

# Apply a taxon filtering to the counts table
process_taxon_selection <- function(ta) {

  # filter counts table
  selected_taxa <- ta$taxa$taxon_id
  ta$counts <- ta$counts %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

# Apply an count filtering to the taxon table
process_count_selection <- function(ta) {

  # filter taxon table
  selected_taxa <- ta$counts$taxon_id %>% unique()
  ta$taxa <- ta$taxa %>%
    filter(taxon_id %in% selected_taxa)

  # return ta object
  ta

}

#' Return rank names associated with a tidytacos object
#'
#' @export
rank_names <- function(ta) {

  if (! is.null(ta$rank_names)) {
    ta$rank_names
  } else {
    c("kingdom", "phylum", "class", "order", "family", "genus")
  }

}

#' Set rank names for a tidytacos object
#'
#' @param ta a tidytacos object
#' @param rank_names a vector containing the names of the ranks used to describe the taxa
#' @return An updated tidytacos object.
#'
#' @export
set_rank_names <- function(ta, rank_names) {

  ta$rank_names <- rank_names

  ta

}

merge_redundant_taxa <- function(ta) {

  # merge taxa in taxon table
  ta$taxa <- ta$taxa %>%
    group_by(taxon_id) %>%
    summarize_all(function(x) {
      x <- unique(x)
      x <- x[! is.na(x)]
      if (length(x) == 1) return(x)
      as.character(NA)
    })

  # merge taxa in counts table
  ta$counts <- ta$counts %>%
    group_by(sample_id, taxon_id) %>%
    summarise(count = sum(count)) %>%
    ungroup()

  ta

}

retain_taxon_id <- function(ta) {
  if ((! "taxon_id" %in% names(ta$taxa)) ||
  (! "taxon_id" %in% names(ta$counts)) ||
  (is.null(ta$counts$taxon_id)) ||
  (is.null(ta$taxa$taxon_id))) {
    stop("You cannot delete the taxon_id column")
  }
}

retain_sample_id <- function(ta) {
  if ((! "sample_id" %in% names(ta$samples)) ||
  (! "sample_id" %in% names(ta$counts)) ||
  (is.null(ta$samples$sample_id)) ||
  (is.null(ta$counts$sample_id))) {
    stop("You cannot delete the sample_id column")
  }
}

retain_counts <- function(ta) {
  if (! "count" %in% names(ta$counts) ||
  (is.null(ta$counts$count))) {
    stop("You cannot delete the count column")
  }
}

any_taxa_left <- function(ta,
  error_message = "No taxa left after filtering") {
    if (length(ta$taxa$taxon_id) == 0) {
    stop(error_message)
  }
}

any_samples_left <- function(ta,
  error_message = "No samples left after filtering") {
    if (length(ta$samples$sample_id) == 0) {
    stop(error_message)
  }
}
