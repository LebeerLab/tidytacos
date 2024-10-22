#' Add relative abundances to count table
#'
#' `add_rel_abundance()` calculates relative abundances of taxa in samples
#' and adds them to the count table of a tidytacos object under the column name
#' "rel_abundance".
#'
#' @param ta A tidytacos object.
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'   c(1500, 1300, 280, 356),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'   taxa_are_columns = FALSE
#' )
#'
#' # Add relative abundance
#' data <- data %>% add_rel_abundance()
#' @family count-modifiers
#' @return A tidytacos object with relative abundances added to the count table.
#' @export
add_rel_abundance <- function(ta) {
  # add relative abundance to abundance table
  ta$counts <- ta$counts %>%
    group_by(sample_id) %>%
    mutate(rel_abundance = count / sum(count)) %>%
    ungroup()

  # return ta object
  ta
}

#' Add absolute abundances to count table
#'
#' `add_absolute_abundance()` calculates absolute abundances of taxa in
#' samples given a taxon that was spiked into all of the samples during library
#' prep. The function then adds these absolute abundances to the count table of
#' the tidytacos object under the column name "absolute_abundance".
#'
#' @param ta A tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name in the sample table which indicates how
#'   much spike was added per sample.
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'   c(1500, 1300, 14, 280, 356, 9),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2", "taxon3")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'   taxa_are_columns = FALSE
#' )
#' data$samples$spike_added <- c(100, 150)
#'
#' # Add total abundance
#' data <- data %>%
#'   add_absolute_abundance(spike_taxon = "t3")
#' @return A tidytacos object with absolute abundances added to the count table.
#' @family count-modifiers
#' @export
add_absolute_abundance <- function(ta, spike_taxon, spike_added = spike_added) {
  spike_added <- rlang::enquo(spike_added)
  if (!rlang::quo_name(spike_added) %in% names(ta$samples)) {
    stop(paste(
      "Sample table requires a column",
      rlang::quo_name(spike_added),
      "that defines the quantity of spike added to the sample."
    ))
  }

  # if total_counts, relative abundances, and total absolute abundances are not
  # present: add temporarily
  total_counts_tmp <- !"total_counts" %in% names(ta$samples)
  if (total_counts_tmp) ta <- add_total_count(ta)
  rel_abundance_tmp <- !"rel_abundance" %in% names(ta$counts)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)
  total_absolute_abundance_tmp <-
    !"total_absolute_abundance" %in% names(ta$samples)
  if (total_absolute_abundance_tmp) ta <-
    add_total_absolute_abundance(ta, spike_taxon)

  # make counts table with total abundances
  tot_abs_abundance <- ta$samples %>%
    select(sample_id, total_absolute_abundance)

  # calculate total for each taxon in each sample
  ta$counts <- ta$counts %>%
    left_join(tot_abs_abundance, by = 'sample_id') %>%
    mutate(
      absolute_abundance =
        round(rel_abundance * total_absolute_abundance, digits = 0)
    )

  # remove total_absolute_abundance from counts
  ta$counts$total_absolute_abundance <- NULL

  # cleanup
  if (total_counts_tmp) ta$samples$total_counts <- NULL
  if (rel_abundance_tmp) ta$counts$rel_abundance <- NULL
  if (total_absolute_abundance_tmp) ta$samples$total_absolute_abundance <- NULL

  # warn about samples without spike
  samples_w_no_spike <-
    unique(ta$counts$sample_id[which(is.na(ta$counts$absolute_abundance))])
  if (length(samples_w_no_spike) > 0) {
    warning(
      paste(
        "Sample without spike taxon detected:",
        format(samples_w_no_spike, trim = TRUE), "\n"
      )
    )
  }

  # return ta object
  ta
}

#' Add density to count table
#'
#' `add_density()` adds densities (bacterial biomass per sample mass or
#' volume) to the count table of a tidytacos object under the column name
#' "density". Can only be used of a taxon was spiked into samples during library prep.
#'
#' @param ta A tidytacos object.
#' @param spike_taxon The taxon id of the spike.
#' @param spike_added The column name of the sample table which indicates how
#'   much spike was added per sample, e.g. 16S rRNA gene copy numbers added to
#'   the DNA extraction tube.
#' @param material_sampled The column name indicating the amount of material
#'   from which DNA was extracted, e.g gram of soil. This parameter encourages
#'   researchers to consider that absolute abundances are only meaningful if
#'   they can be translated into densities.
#'
#' @examples
#' # Initiate count matrix
#' x <- matrix(
#'   c(1500, 1300, 14, 280, 356, 9),
#'   ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2", "taxon3")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'   taxa_are_columns = FALSE
#' )
#' data$samples$spike_added <- c(100, 150)
#' data$samples$grams_source <- c(4, 5)
#'
#' # Add density
#' data <- data %>%
#'   add_density(spike_taxon = "t3", material_sampled=grams_source)
#' @return A tidytacos object with densities added to the count table.
#' @family count-modifiers
#' @export
add_density <- function(
    ta, spike_taxon, spike_added = spike_added,
    material_sampled = material_sampled
  ) {

  spike_added <- rlang::enquo(spike_added)
  material_sampled <- rlang::enquo(material_sampled)

  if (!rlang::as_name(spike_added) %in% names(ta$samples)) {
    stop(paste(
      "Sample table requires a column",
      rlang::as_name(spike_added),
      "that defines the quantity of spike added to the sample."
    ))
  }

  if (!rlang::as_name(material_sampled) %in% names(ta$samples)) {
    stop(paste(
      "Sample table requires a column",
      rlang::as_name(material_sampled),
      "that defines the quantity of sample used."
    ))
  }

  # if relative abundances and total densities are not present: add temporarily
  rel_abundance_tmp <- !"rel_abundance" %in% names(ta$counts)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)
  total_density_tmp <- !"total_density" %in% names(ta$samples)
  if (total_density_tmp) ta <- add_total_density(ta, spike_taxon,
          spike_added = !!spike_added, material_sampled = !!material_sampled)

  # make counts table with total densities
  tot_densities <- ta$samples %>%
    select(sample_id, total_density)

  # calculate density for each taxon in each sample
  ta$counts <- ta$counts %>%
    left_join(tot_densities, by='sample_id') %>%
    mutate(density = round( rel_abundance * total_density, digits = 0))

  # remove total_density from counts
  ta$counts$total_density <- NULL

  # cleanup
  if (rel_abundance_tmp) ta$counts$rel_abundance <- NULL
  if (total_density_tmp) ta$samples$total_density <- NULL

  # Warn about samples without spike
  samples_w_no_spike <-
    unique(ta$counts$sample_id[which(is.na(ta$counts$density))])
  if (length(samples_w_no_spike) > 0) {
    warning(
      paste(
        "Sample without spike taxon detected:",
        format(samples_w_no_spike, trim = TRUE), "\n"
      )
    )
  }

  # return ta object
  ta
}
