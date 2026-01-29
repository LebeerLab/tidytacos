#' Rarefy the samples to a given number of reads
#'
#' This function performs rarefying. Make sure that all samples contain at least
#' the minimum number of reads; otherwise, an error might be thrown.
#'
#' @param ta A tidytacos object.
#' @param n Subsample size for rarefying the community.
#' @param replace Whether to replace the read after it has been selected for the subsample so it can be sampled again. Default is FALSE.
#' @return A tidytacos object.
#'
#' @examples
#' # discard samples with less than 1000 reads
#' urt_1000 <- urt %>%
#'   add_total_count() %>%
#'   filter_samples(total_count >= 1000)
#'
#' # then rarefy to 1000 reads
#' urt_1000 <- urt_1000 %>% rarefy(1000)
#'
#' @export
rarefy <- function(ta, n, replace = F) {
  ta$counts <- try(
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
  )

  if (class(ta$counts)[[1]] == "try-error") {
    stop(
      paste(
        "Rarefying failed. Make sure that all samples contain at least the minimum number of reads.\n",
        "Or use replace = TRUE, to allow sampling with replacement."
      )
    )
  }

  ta %>%
    purrr::modify_at("counts", filter, count > 0) %>%
    process_count_selection()
}

#' Change sample IDs to a given expression
#'
#' @param ta A tidytacos object.
#' @param sample_id_new An expression that evaluates to a unique sample
#'   identifier.
#'
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

#' Change taxon IDs to a given expression
#'
#' @param ta A tidytacos object.
#' @param taxon_id_new An expression that evaluates to a unique taxon
#'   identifier.
#'
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
#' `aggregate_samples()` merges sample content of samples which have identical values for all columns in the sample table (except sample_id).
#'
#' @param ta A tidytacos object.
#'
#' @return A tidytacos object.
#' @export
aggregate_samples <- function(ta) {
  sample_id_new <- NULL

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
#' @param ta A tidytacos object.
#' @param rank An optional rank to aggregate on.
#' @return A tidytacos object.
#'
#' @examples
#' urt %>% aggregate_taxa(rank = "class")
#' @export
aggregate_taxa <- function(ta, rank = NULL) {
  taxon_id_new <- NULL

  # Temporarily replace any NA's with strings as they interfere with aggregation
  ta$taxa[is.na(ta$taxa)] <- "unknown"

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
    ta <- suppressWarnings(select_taxa(ta, taxon_id, !!rank_names_to_keep))
  }

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
  # Adapt rank names to aggregate
  ta <- ta %>% set_rank_names(
    rank_names(ta) %>% intersect(names(ta$taxa))
  )
  # Add new unique taxon label
  if (!is.null(rank)) {
    include_species <- eval(rank == "species")
    ta <- ta %>%
      add_taxon_name(include_species = include_species) %>%
      mutate_taxa(taxon = taxon_name) %>%
      suppressWarnings(select_taxa(-taxon_name))
  }
  ta
}

#' Trim all sequences
#'
#' `trim_asvs()` trims sequence ends of the sequence supplied in the taxa table.
#' This function assumes that the sequence variable in the taxon table is called
#' "sequence".
#'
#' @param ta A tidytacos object.
#' @param start Index of where to start trimming.
#' @param end Index of where to stop trimming.
#'
#' @return A tidytacos object.
#' @examples
#' # keep only the first 200 nucleotides of the sequences
#' urt %>% trim_asvs(0, 200)
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

  ta <- aggregate_trimmed_asvs(ta)
  ta <- merge_redundant_taxa(ta)

  ta
}

#' Retain or remove a set of sample variables
#'
#' @param ta A tidytacos object.
#' @param ... Selection criteria for the samples table.
#' @return A tidytacos object.
#'
#' @examples
#' # remove the condition column from the samples table
#' urt %>% select_samples(-condition)
#' # keep only the sample_id, location and method columns
#' urt %>% select_samples(sample_id, location, method)
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
#'
#' @param ta A tidytacos object.
#' @param ... Selection criteria for the taxa table.
#' @return A tidytacos object.
#' @examples
#'
#' # drop the sequence column
#' urt %>% select_taxa(-sequence)
#'
#' # keep only the taxon_id and genus columns
#' urt %>% select_taxa(taxon_id, genus)
#'
#' @export
select_taxa <- function(ta, ...) {
  ta$taxa <- ta$taxa %>%
    select(...)

  rn_missing <- setdiff(ta %>% rank_names(), colnames(ta$taxa))
  rn <- ta %>% rank_names()

  retain_taxon_id(ta)

  ta %>% set_rank_names(rn[!rn %in% rn_missing])
}

#' Retain or remove a set of count variables
#'
#' @param ta A tidytacos object.
#' @param ... Selection criteria for the counts table.
#' @return A tidytacos object.
#' @examples
#'
#' # add a column to the counts table
#' leaf_ab <- leaf %>% add_rel_abundance()
#' # remove that column again
#' leaf_ab %>% select_counts(-rel_abundance)
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
#'
#' @param ta A tidytacos object.
#' @param ... Mutate criteria for the samples table.
#' @return A tidytacos object.
#' @examples
#'
#' # change the sample column to lowercase
#' urt <- urt %>% mutate_samples(sample = tolower(sample))
#'
#' @export
mutate_samples <- function(ta, ...) {
  ta$samples <- ta$samples %>%
    mutate(...)
  retain_sample_id(ta)

  ta
}

#' Create extra variables in the taxa table
#'
#' @param ta A tidytacos object.
#' @param ... Mutate criteria for the taxa table.
#' @return A tidytacos object.
#' @examples
#' urt <- urt %>% mutate_taxa(species = paste(genus, species))
#' @export
mutate_taxa <- function(ta, ...) {
  ta$taxa <- ta$taxa %>%
    mutate(...)
  retain_taxon_id(ta)

  ta
}

#' Create extra variables in the count table
#'
#' @param ta A tidytacos object.
#' @param ... Mutate criteria for the counts table.
#' @return A tidytacos object.
#' @examples
#' # add a column to the counts table
#' urt %>% mutate_counts(cl10 = log10(count))
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
#'
#' @param ta A tidytacos object.
#' @param ... Filter criteria for the samples table.
#' @return A tidytacos object.
#' @examples
#'
#' # subset urt to keep only nasopharynx samples
#' urt_nf <- urt %>% filter_samples(location == "NF")
#' # subset urt to keep only samples from plate 1 and 2
#' urt_plate_1_2 <- urt %>% filter_samples(plate %in% c(1, 2))
#' # subset the blanks in leaf
#' leaf_blanks <- leaf %>% filter_samples(startsWith(description, "BLANK"))
#'
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
#'
#' @param ta A tidytacos object.
#' @param ... Filter criteria for the taxa table.
#' @return A tidytacos object.
#' @examples
#' # keep only bacterial reads
#' leaf <- leaf %>% filter_taxa(kingdom == "Bacteria")
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
#'
#' @param ta A tidytacos object.
#' @param ... Filter criteria for the counts table.
#' @return A tidytacos object.
#' @examples
#' # remove singletons
#' urt <- urt %>% filter_counts(count > 1)
#' @export
filter_counts <- function(ta, ...) {
  ta$counts <- ta$counts %>%
    filter(...)

  ta <- ta %>%
    process_count_selection()
  any_taxa_left(ta)

  ta
}


#' Group the samples
#'
#' @param ta A tidytacos object.
#' @param ... Grouping criteria for the samples table.
#' @return A (named) list of tidytacos object.
#' @examples
#'
#' urt_by_loc <- urt %>% group_samples(location)
#' 
#' # apply a function to each separate taco, eg. tacosum
#' urt_by_loc@tacos %>% lapply(tacosum)
#' 
#' # subset urt to keep only nasopharynx samples
#' urt_nf <- urt_by_loc@tacos$NF
#' @export
group_samples <- function(ta, ...) {
    gr_names <- ta$samples %>%
      dplyr::group_by(...) %>%
      dplyr::group_keys()
  
    gr_cols <- colnames(gr_names)

    sample_groups <- ta$samples %>%
      dplyr::group_split(...)
    groups <- lapply(sample_groups, function (x) ta %>% filter_samples(sample_id %in% x$sample_id))

    if (length(gr_cols) == 1) {
        names(groups) <- gr_names %>% dplyr::pull()
    } else {
      names(groups) <- gr_names %>% tidyr::unite(label) %>% pull(label)
    }
    new("grouped_taco", tacos=groups, group_cols=gr_cols)
}

#' Perform a centered log ratio transformation on the readcounts.
#'
#' `add_clr_abundance()` calculates the log ration transformed values for each taxon in each sample and adds these data in a new table, clr_counts. Alternatively, using 'overwrite', the clr transformed data can replace the 'counts' column in the count table.
#'
#' @param ta A tidytacos object.
#' @param overwrite Whether or not the counts table is to be overwritten with the transformed counts.
#' @param pseudocount A pseudocount to be added to the counts before transformation. If false or zero will perform robust CLR.
#' @inheritDotParams vegan::decostand
#' @return A tidytacos object.
#' @export
add_clr_abundance <- function(
    ta,
    overwrite = F,
    pseudocount = 1,
    ...) {
  
  mt <- ta %>% counts_matrix()

  if (pseudocount == 0) {
    clrt_mt <- vegan::decostand(mt, method="rclr", ...)
  } else {
    clrt_mt <- vegan::decostand(mt, method="clr", pseudocount=pseudocount, ...)
  }

  clrtt <- clrt_mt %>% create_tidytacos(allow_non_count=TRUE)
  counts_w_prev_indexers <- clrtt %>% 
    everything() %>%
    select(sample, taxon, count) %>%
    rename(taxon_id=taxon, sample_id = sample)   
 
  if (overwrite) {
    ta$counts <- counts_w_prev_indexers
  } else {
    if (pseudocount == 0) {
      ta$rclr_counts <- counts_w_prev_indexers
    } else {
      ta$clr_counts <- counts_w_prev_indexers
    }
  }

  ta
}
