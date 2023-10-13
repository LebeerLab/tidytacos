#' Rarefy the samples to a given number of reads
#'
#' This function performs rarefying. Make sure that all samples contain at least
#' the minimum number of reads; otherwise, an error might be thrown.
#'
#' @export
rarefy <- function(ta, n, replace = F) {
  ta$counts <-
    ta$counts %>%
    group_by(sample_id) %>%
    mutate(
      count =
        sample(x = 1:sum(count), size = !!n, replace = !!replace) %>%
          cut(breaks = c(0, cumsum(count)), labels = taxon_id) %>%
          table() %>%
          as.integer()
    ) %>%
    ungroup()

  ta %>%
    purrr::modify_at("counts", filter, count > 0) %>%
    process_count_selection()
}

# Change sample IDs to a given expression
#
# @param ta A tidytacos object.
# @param sample_id_new An expression that evaluates to a unique sample
#   identifier.
#
change_id_samples <- function(ta, sample_id_new) {
  sample_id_new <- rlang::enexpr(sample_id_new)

  ta <- mutate_samples(ta, sample_id_new = as.character(!!sample_id_new))

  if (any(duplicated(ta$samples$sample_id_new))) {
    stop("the new sample ids are not unique")
  }

  ta$counts <-
    ta$counts %>%
    left_join(
      ta$samples %>% select(sample_id, sample_id_new),
      by = "sample_id"
    ) %>%
    select(-sample_id) %>%
    rename(sample_id = sample_id_new)

  ta$samples <-
    ta$samples %>%
    select(-sample_id) %>%
    rename(sample_id = sample_id_new)

  ta
}

# Change taxon IDs to a given expression
#
# @param ta A tidytacos object.
# @param taxon_id_new An expression that evaluates to a unique taxon
#   identifier.
#
change_id_taxa <- function(ta, taxon_id_new) {
  taxon_id_new <- rlang::enexpr(taxon_id_new)

  ta <- mutate_taxa(ta, taxon_id_new = as.character(!!taxon_id_new))

  if (any(duplicated(ta$taxa$taxon_id_new))) {
    stop("the new taxon ids are not unique")
  }

  ta$counts <-
    ta$counts %>%
    left_join(ta$taxa %>% select(taxon_id, taxon_id_new), by = "taxon_id") %>%
    select(-taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta$taxa <-
    ta$taxa %>%
    select(-taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta
}

#' Aggregate samples with identical values for all metadata
#'
#' @param ta A tidytacos object.
#'
#' @export
aggregate_samples <- function(ta) {
  # sample table with only old and new sample names
  metadata <- setdiff(names(ta$samples), "sample_id")
  names <- ta$samples %>%
    select(-sample_id) %>%
    distinct() %>%
    mutate(sample_id_new = paste0("m", 1:n())) %>%
    right_join(ta$samples, by = metadata, multiple = "all") %>%
    select(sample_id, sample_id_new)

  # adapt sample table with new names
  ta$samples <- ta$samples %>%
    left_join(names, by = "sample_id") %>%
    select(-sample_id) %>%
    rename(sample_id = sample_id_new) %>%
    distinct()

  # merge samples in counts table and adapt with new names
  ta$counts <- ta$counts %>%
    left_join(names, by = "sample_id") %>%
    select(-sample_id) %>%
    group_by(sample_id_new, taxon_id) %>%
    summarize(count = sum(count)) %>%
    ungroup() %>%
    rename(sample_id = sample_id_new)

  # return ta object
  ta
}

#' Aggregate taxa on a given taxonomic rank
#'
#' There are two ways to call this function:
#'
#' * If the rank you are interested in is in the standard list, just supply it
#' as an argument.
#' * If not, delete all taxon variables except taxon_id and the ranks you are
#' still interested in prior to calling this function.
#'
#' @param ta a tidytacos object
#' @param rank an optional rank to aggregate on
#' @export
aggregate_taxa <- function(ta, rank = NULL) {
  if (!is.null(rank)) {
    rank_names <-
      rank_names(ta) %>%
      intersect(names(ta$taxa))

    if (length(rank_names) == 0) {
      stop(
        "at least one of the taxonomic rank names should be present ",
        "in the taxon table"
      )
    }

    if (!rank %in% rank_names) {
      stop("the rank you supplied should be one of the rank names")
    }

    rank_index <- which(rank_names == rank)
    rank_names_to_keep <- rank_names[1:rank_index]
    ta <- select_taxa(ta, taxon_id, !!rank_names_to_keep)
  }

  # this avoids some problems
  ta$taxa[is.na(ta$taxa)] <- "unknown"

  ta$taxa <-
    ta$taxa %>%
    chop(taxon_id) %>%
    mutate(taxon_id_new = paste0("t", 1:n()))

  id_conversion <-
    ta$taxa %>%
    unnest(taxon_id) %>%
    select(taxon_id, taxon_id_new)

  ta$taxa <-
    ta$taxa %>%
    select(-taxon_id) %>%
    rename(taxon_id = taxon_id_new)

  ta$counts <-
    ta$counts %>%
    left_join(id_conversion, by = "taxon_id") %>%
    select(-taxon_id) %>%
    group_by(taxon_id_new, sample_id) %>%
    {
      if ("rel_abundance" %in% names(ta$counts)) {
        summarize(
          .,
          count = sum(count), rel_abundance = sum(rel_abundance)
        )
      } else {
        summarize(., count = sum(count))
      }
    } %>%
    ungroup() %>%
    rename(taxon_id = taxon_id_new)

  # cleanup
  ta$taxa[ta$taxa == "unknown"] <- NA

  ta
}

#' Trim all sequences
#'
#' This function assumes that the sequence variable in the taxon table is called
#' "sequence".
#' @param ta a tidytacos object
#' @param start index of where to start trimming
#' @param end index of where to stop trimming
#'
#' @export
trim_asvs <- function(ta, start, end) {
  ta$taxa <- ta$taxa %>%
    mutate(sequence = str_sub(sequence, start = !!start, end = !!end))
  if ("sequence" %in% names(ta$counts)) {
    ta$counts <- ta$counts %>%
      mutate(sequence = str_sub(
        sequence,
        start = !!start, end = !!end
      ))
  }
  ta <- merge_redundant_taxa(ta)

  ta
}

#' Retain or remove a set of sample variables
#' @param ta a tidytacos object
#' @export
select_samples <- function(ta, ...) {
  ta$samples <- ta$samples %>%
    select(...)

  if (!"sample_id" %in% names(ta$samples)) {
    stop("you cannot delete the sample_id column")
  }

  ta
}

#' Retain or remove a set of taxon variables
#' @param ta a tidytacos object
#' @export
select_taxa <- function(ta, ...) {
  ta$taxa <- ta$taxa %>%
    select(...)

  retain_taxon_id(ta)

  ta
}

#' Retain or remove a set of count variables
#' @param ta a tidytacos object
#' @export
select_counts <- function(ta, ...) {
  ta$counts <- ta$counts %>%
    select(...)

  retain_sample_id(ta)
  retain_taxon_id(ta)
  retain_counts(ta)

  ta
}

#' Create extra variables in the sample table
#' @param ta a tidytacos object
#' @export
mutate_samples <- function(ta, ...) {
  ta$samples <- ta$samples %>%
    mutate(...)
  retain_sample_id(ta)

  ta
}

#' Create extra variables in the taxon table
#' @param ta a tidytacos object
#' @export
mutate_taxa <- function(ta, ...) {
  ta$taxa <- ta$taxa %>%
    mutate(...)
  retain_taxon_id(ta)

  ta
}

#' Create extra variables in the abundances table
#' @param ta a tidytacos object
#' @export
mutate_counts <- function(ta, ...) {
  ta$counts <- ta$counts %>%
    mutate(...)
  retain_sample_id(ta)
  retain_taxon_id(ta)
  retain_counts(ta)

  ta
}

#' Filter the samples
#' @param ta a tidytacos object
#' @export
filter_samples <- function(ta, ...) {
  ta$samples <- ta$samples %>%
    filter(...)

  ta <- ta %>%
    process_sample_selection()
  any_samples_left(ta)

  ta
}

#' Filter the taxa
#' @param ta a tidytacos object
#' @export
filter_taxa <- function(ta, ...) {
  ta$taxa <- ta$taxa %>%
    filter(...)

  ta <- ta %>%
    process_taxon_selection()
  any_taxa_left(ta)

  ta
}

#' Filter the counts
#' @param ta a tidytacos object
#' @export
filter_counts <- function(ta, ...) {
  ta$counts <- ta$counts %>%
    filter(...)

  ta <- ta %>%
    process_count_selection()
  any_taxa_left(ta)

  ta
}

#' Perform a centered log ratio transformation on the readcounts.
#' @param ta a tidytacos object
#' @param taxon unique identifier for the taxon, by default taxon_id.
#' @param sample unique identifier for the sample, by default sample_id.
#' @param counts variable name of the counts to be transformed in the counts table.
#' @param overwrite wether or not the counts table is to be overwritten with the transformed counts.
#' @export
clr_transform_counts <- function(
    ta,
    overwrite = F) {
  force_optional_dependency("compositions")

  mt <- ta$counts %>% pivot_wider(
    values_from=count,
    names_from=taxon_id,
    values_fill=0)

  mt <- tibble::column_to_rownames(mt, var="sample_id")

  clrt_mt <- compositions::clr(mt) 
  clrt_counts <- clrt_mt %>% 
      as_tibble() %>% 
      tibble::add_column(sample_id = rownames(clrt_mt)) %>% 
      pivot_longer(!sample_id, names_to='taxon_id', values_to='count') %>% 
      filter(count>0)

  if (overwrite) {
    ta$counts <- clrt_counts
  } else {
    ta$clr_counts <- clrt_counts
  }

  ta
}
