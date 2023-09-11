#' (Re)classify amplicon sequences
#'
#' This function requires the DADA2 package to be installed.
#'
#' This function will (re)classify either all or a subset of the taxa, given
#' that a variable is present in the taxon table that contains (representative)
#' sequences of the taxa.
#'
#' Ranks can be supplied as a named integer vector, where the names represent
#' taxonomic ranks and the integers represent positions of these ranks in the
#' taxonomy strings present in the reference database. Ranks can also be
#' supplied as just a character vector with the rank names; in that case, it is
#' assumed that the database taxonomy string follows the default order {domain,
#' phylum, class, order, family, genus, species}. If no ranks are supplied, taxa
#' will be (re)classified at all default ranks.
#'
#' @param ta A tidytacos object.
#' @param refdb The path to a DADA2-compatible reference database.
#' @param taxa An expression that specifies which taxa to (re)classify.
#' @param ranks A vector that specifies which ranks to (re)classify.
#' @param sequence_var The (quoted) name of a variable within the taxa table
#'   that contains (representative) sequences of the taxa.
#' @param multithread A boolean indicating whether to use multiple threads.
#' @param min_boot The minimum bootstrap value for taxonomy assignment.
#' @param n_ranks The number of ranks present in the reference database.
#'
#' @return An updated tidytacos object.
#'
#' @export
classify_taxa <- function(
  ta, refdb, taxa = rep(T, times = length(taxon_id)), ranks = "default",
  sequence_var = "sequence", multithread = T, min_boot = 50, n_ranks = 7
) {

  force_optional_dependency("dada2")
  # throw error if sequence_var doesn't exist
  if (! sequence_var %in% names(ta$taxa)) {
    stop(paste0("variable '", sequence_var, "' not found in taxon table"))
  }

  # convert taxa argument to an expression
  taxa <- rlang::enexpr(taxa)

  # determine which taxa to (re)classify
  to_reclassify <- eval(taxa, ta$taxa)
  to_reclassify <- to_reclassify & ! is.na(to_reclassify)

  # throw error if no taxa to (re)classify
  if (sum(to_reclassify) == 0) stop("zero taxa obey the given criteria")

  # lookup default ranks and/or rank numbers if not given
  if (! is.numeric(ranks)) {
    default_ranks <- 1:7
    names(default_ranks) <-
      c("domain", "phylum", "class", "order", "family", "genus", "species")
    if (ranks[1] == "default") ranks <- names(default_ranks)
    ranks <- default_ranks[ranks]
    ranks <- ranks[! is.na(ranks)]
  }

  # throw error if a rank number exceeds the number of ranks in the db
  if (max(ranks) > n_ranks) stop("not enough ranks in the database")

  # perform (re)classification
  # tryRC is important; see <https://github.com/benjjneb/dada2/issues/1441>
  taxa_reclassified <-
    ta$taxa[to_reclassify, ][[sequence_var]] %>%
    dada2::assignTaxonomy(
      refFasta = refdb, multithread = multithread, minBoot = min_boot,
      tryRC = T, taxLevels = letters[1:n_ranks]
    )
  ta$taxa[to_reclassify, names(ranks)] <- taxa_reclassified[ , ranks]

  ta

}

#' Add taxon metadata to the tidytacos object
#'
#' \code{add_taxon_tibble} adds a taxon tibble to the tidytacos object.
#'
#' This function adds a taxon tibble containing taxon data (e.g. taxon ranks
#' such as genus, family, ...) for each taxon to the tidytacos object. It is
#' used after initiating a tidytacos object using a numerical abundance
#' matrix and the function \code{\link{tidytacos}}. Also see
#' \code{\link{add_sample_tibble}} to update the sample data of the
#' tidytacos object.
#'
#' @param ta tidytacos object.
#' @param taxon_tibble A tibble containing taxon data for each taxon. Taxa
#'   should be rows, while taxon data should be columns. At least one column
#'   name needs to be shared with the taxon tibble of ta. The default shared
#'   column name is 'taxon'.
#'
#' @examples
#' # Initiate abundance matrix
#' x <- matrix(
#'  c(1500, 1300, 280, 356),
#'  ncol = 2
#' )
#' rownames(x) <- c("taxon1", "taxon2")
#' colnames(x) <- c("sample1", "sample2")
#'
#' # Convert to tidytacos object
#' data <- create_tidytacos(x,
#'                      taxa_are_columns = FALSE)
#'
#' # Initiate taxon tibble
#' taxon <- c("taxon1", "taxon2")
#' genus <- c("Salmonella", "Lactobacillus")
#' taxon_tibble <- tibble::tibble(taxon, genus)
#'
#' # Add taxon tibble to tidytacos object
#' data <- data %>%
#' add_taxon_tibble(taxon_tibble)
#'
#' @export
add_taxon_tibble <- function(ta, taxon_tibble) {

  purrr::modify_at(ta, "taxa", left_join, taxon_tibble)

}

#' Create sensible names for the taxa and add to taxon table
#'
#' @param ta A tidytacos object.
#' @param method The method by which to arrange the taxon names. Currently only
#'   mean_rel_abundance.
#' @param include_species Whether to include the species name or not.
#'
#' @importFrom stats na.omit
#' @export
add_taxon_name <- function(
  ta, method = "mean_rel_abundance", include_species = F
  ) {

  if (method == "mean_rel_abundance") {

    # if mean_rel_abundance not present: add temporarily
    mean_rel_ab_tmp <- ! "mean_rel_abundance" %in% names(ta$taxa)
    if (mean_rel_ab_tmp) ta <- add_mean_rel_abundances(ta)

    ta <- mutate_taxa(ta, arrange_by_me = mean_rel_abundance)

  } else {

    # throw error if method unknown
    stop("method unknown")

  }

  rank_names <- rank_names(ta)
  if (include_species) {
    rank_names <- append(rank_names, "species", after=length(rank_names))
  }

  if (length(rank_names) == 0) {

    ta <- mutate_taxa(ta, best_classification = "unclassified")

  } else {

    ta$taxa <-
      ta$taxa %>%
      mutate(
        best_classification =
          purrr::pmap_chr(
            ta$taxa[, rank_names],
            function(...) {
              classification = as.character(list(...))
              if (all(is.na(classification))) return("unclassified")
              classification %>% na.omit() %>% last()
            }
          )
      )

  }

  ta$taxa <-
    ta$taxa %>%
    group_by(best_classification) %>%
    arrange(desc(arrange_by_me)) %>%
    mutate(n_taxa = n()) %>%
    mutate(taxon_number = if_else(n_taxa > 1, as.character(1:n()), "")) %>%
    mutate(taxon_name = str_c(best_classification, taxon_number, sep = " ")) %>%
    mutate_at("taxon_name", str_trim) %>%
    ungroup() %>%
    select(- best_classification, - n_taxa, - taxon_number)

  # cleanup
  if (exists("mean_rel_ab_tmp")) ta$taxa$mean_rel_abundance <- NULL
  ta$taxa$arrange_by_me <- NULL

  # return ta object
  ta

}

#' Create taxon names suitable for visualization with color.
#'
#' A rank can be supplied to aggregate colors higher than the current rank.
#'
#' @param ta A tidytacos object.
#' @param method The method by which to arrange the taxon names. Currently only
#'   mean_rel_abundance.
#' @param n An integer denoting the amount of most abundant taxa to display.
#'   Capacity at 12.
#' @param samples An optional vector of sample_id's of interest.
#' @param taxa An optional vector of taxon_id's of interest.
#' @param rank An optional rank to aggregate taxa on.
#'
#' @export
add_taxon_name_color <- function(
  ta, method = "mean_rel_abundance", n = 12, samples = NULL, taxa = NULL,
  rank = NULL
  ) {

  # Recursive function to add taxon name color on different rank
  add_colnames_rank <- function(ta, rank_col) {
    colnames <- ta %>%
      aggregate_taxa(rank = rank_col) %>%
      add_taxon_name_color(method=method, n=n)

    col_per_sample <- colnames$counts %>%
      inner_join(colnames$taxa, by="taxon_id") %>%
      dplyr::select(sample_id, taxon_name_color, taxon_id)

    ta$taxa <- ta$taxa %>% left_join(col_per_sample, by="taxon_id")
    ta

  }

  # if taxon_name not present: add temporarily
  taxon_name_tmp <- ! "taxon_name" %in% names(ta$taxa)
  if (taxon_name_tmp) ta <- add_taxon_name(ta)

  if (method == "mean_rel_abundance") {

    # if mean_rel_abundance not present: add temporarily
    mean_rel_ab_tmp <- ! "mean_rel_abundance" %in% names(ta$taxa)
    if (mean_rel_ab_tmp) ta <- add_mean_rel_abundances(ta)

    ta <- mutate_taxa(ta, arrange_by_me = mean_rel_abundance)

  } else {

    # throw error if method unknown
    stop("method unknown")

  }

  ta_subset <- ta

  # take subset of samples if requested
  if (! is.null(samples)) {
    ta_subset <- filter_samples(ta_subset, sample_id %in% !! samples)
  }

  # take subset of taxa if requested
  if (! is.null(taxa)) {
    ta_subset <- filter_taxa(ta_subset, taxon_id %in% !! taxa)
  }

  # extract taxon names to visualize, in order
  levels <-
    ta_subset$taxa %>%
    arrange(desc(arrange_by_me)) %>%
    pull(taxon_name) %>%
    `[`(1:(n-1)) %>%
    sort()
  levels <- append(levels, "residual", after=0)

  # add taxon_name_color factor to taxa table
  ta$taxa <-
    ta$taxa %>%
    mutate(taxon_name_color = if_else(taxon_name %in% levels, taxon_name, "residual")) %>%
    mutate(taxon_name_color = factor(taxon_name_color, levels = levels))

  # overwrite colors with higher rank if asked for
  if(! is.null(rank)) {
    ta$taxa <- ta$taxa %>% dplyr::select(-one_of("taxon_name_color"))
    ta <- add_colnames_rank(ta, rank)
  }

  # cleanup
  if (taxon_name_tmp) ta$taxa$taxon_name <- NULL
  if (exists("mean_rel_ab_tmp")) ta$taxa$mean_rel_abundance <- NULL
  ta$taxa$arrange_by_me <- NULL

  # return ta object
  ta

}

#' Apply the taxon QC method of Jervis-Bardy
#'
#' Function to estimate spearman correlation between relative abundance and
#' sample dna concentration, for each taxon.
#'
#' See:
#' J. Jervis-Bardy et al., “Deriving accurate microbiota profiles from
#' human samples with low bacterial content through post-sequencing processing
#' of Illumina MiSeq data,” Microbiome, vol. 3, no. 1, Art. no. 1, 2015, doi:
#' 10.1186/s40168-015-0083-8.
#'
#' @importFrom stats cor.test
#' @param ta A tidytacos object.
#' @param dna_conc A variable in the samples table that contains dna
#'   concetrations (unquoted).
#' @param sample_condition An optional extra condition that samples must pass
#'   before calculations.
#' @param min_pres The minimum number of samples a taxon has to be present in
#'   for its correlation to be calculated.
#'
#' @export
add_jervis_bardy <- function(ta, dna_conc, sample_condition = T, min_pres = 3) {

  dna_conc <- enquo(dna_conc)
  sample_condition <- enquo(sample_condition)

  # if rel_abundance not present: add temporarily
  rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$counts)
  if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

  # if sample condition is given, use only samples that fulfill it
  if (is.null(sample_condition)) {
    ta_jb <- ta
  } else {
    ta_jb <- ta %>%
      filter_samples(!! sample_condition)
  }

  # perform jervis bardy calculation
  taxa_jb <- ta_jb$counts %>%
    left_join(
      ta_jb$samples %>% select(sample_id, dna_conc = !! dna_conc),
      by = "sample_id"
    ) %>%
    group_by(taxon_id) %>%
    filter(n() >= !! min_pres) %>%
    do(
      jb = cor.test(
        x = .$rel_abundance, y = .$dna_conc, alternative = "less",
        method = "spearman"
      )
    ) %>%
    mutate(jb_cor = jb$estimate, jb_p = jb$p.value) %>%
    select(- jb)

  # add jb_p and jb_cor to taxa table
  ta$taxa <- left_join(ta$taxa, taxa_jb, by = "taxon_id")

  # cleanup
  if (rel_abundance_tmp) ta$counts$rel_abundance <- NULL

  # return ta object
  ta

}

#' Add taxon occurrences to the taxon table
#'
#' Adds taxon occurrences (overall or per condition) to the taxa table.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#' @importFrom stats fisher.test
#' @param ta A tidytacos object
#' @param condition A categorical variable (string)
#' @param relative wether to use relative occurences
#' @param fischer_test wether to perform a fischer test and add the p-values of the test to the taxa table
#' @export
add_occurrences <- function(
  ta, condition = NULL, relative = F, fischer_test = F
  ) {

  if (is.null(condition)) {

    taxa_occurrences <-
      occurrences(ta, condition = condition)

  } else if (fischer_test) {

    occurrences <-
      occurrences(ta, condition = condition, pres_abs = T)

    condition_sym <- sym(condition)

    taxa_fischer <-
      occurrences %>%
      group_by(taxon_id) %>%
      arrange(!! condition_sym, presence) %>%
      do(
        fisher = c(.$n) %>%
          matrix(ncol = 2, byrow = T) %>%
          fisher.test()
      ) %>%
      mutate(fisher_p = fisher$p.value) %>%
      select(- fisher)

    taxa_occurrences <-
      occurrences %>%
      filter(presence == "present") %>%
      select(taxon_id, !! condition_sym, occurrence = n) %>%
      mutate_at(condition, ~ str_c("occurrence_in", ., sep = "_")) %>%
      spread(value = "occurrence", key = condition) %>%
      left_join(taxa_fischer, by = "taxon_id")

  } else {

    taxa_occurrences <-
      occurrences(ta, condition = condition) %>%
      mutate_at(condition, ~ str_c("occurrence_in", ., sep = "_")) %>%
      spread(value = "occurrence", key = condition)

  }

  if (relative & is.null(condition)) {

    taxa_occurrences <-
      taxa_occurrences %>%
      mutate(occurrence = occurrence / nrow(ta$samples))

  }

  if (relative & ! is.null(condition)) {

    condition_sym <- sym(condition)

    conditions <-
      ta$samples %>%
      count(!! condition_sym)

    for (con_ix in 1:nrow(conditions)) {

      con <- conditions[[condition]][con_ix]
      n_samples <- conditions[["n"]][con_ix]
      taxa_occurrences[[str_c("occurrence_in_", con)]] <-
        taxa_occurrences[[str_c("occurrence_in_", con)]] / n_samples

    }

  }

  ta %>%
    purrr::modify_at("taxa", left_join, taxa_occurrences, by = "taxon_id")

}

#' Add average relative abundances
#'
#' This function adds mean relative abundance values for each taxon to the taxa
#' table, overall or per sample group.
#'
#' If `condition` is specified, the mean relative abundances will be calculated
#' separately for each group defined by the condition variable. This variable
#' should be present in the sample table.
#'
#' If `condition` is specified, differential abundance testing can be performed
#' by setting the `test` argument. Options are NULL (default), "wilcox" or
#' "t-test".
#'
#' @importFrom stats t.test wilcox.test
#' @param ta A tidytacos object
#' @param condition A condition variable (character)
#' @param test Differential abundance test to perform
#'
#' @return A tidytacos object
#'
#' @export
add_mean_rel_abundances <- function(ta, condition = NULL, test = NULL) {

  mean_rel_abundances <- mean_rel_abundances(ta, condition = condition)

  if (is.null(condition)) {

    taxa_mean_rel_abundances <- mean_rel_abundances

  } else if (! is.null(test)) {

    condition_sym <- ensym(condition)

    # if rel_abundance not present: add temporarily
    rel_abundance_tmp <- ! "rel_abundance" %in% names(ta$counts)
    if (rel_abundance_tmp) ta <- add_rel_abundance(ta)

    rel_abundances_complete <-
      counts(ta) %>%
      left_join(ta$samples, by = "sample_id") %>%
      select(sample_id, !! condition_sym, taxon_id, rel_abundance) %>%
      complete(
        nesting(sample_id, !! condition_sym), taxon_id,
        fill = list(rel_abundance = 0)
      )

    if (test == "wilcox") {
      taxa_test <-
        rel_abundances_complete %>%
        group_by(taxon_id) %>%
        do(result = wilcox.test(
          data = ., rel_abundance ~ !! condition_sym
        )) %>%
        mutate(wilcox_p = result$p.value, wilcox_stat = result$statistic) %>%
        select(- result)
    } else if (test == "t-test") {
      taxa_test <-
        rel_abundances_complete %>%
        group_by(taxon_id) %>%
        do(result = t.test(
          data = ., rel_abundance ~ !! condition_sym
        )) %>%
        mutate(t_test_p = result$p.value, t_test_stat = result$statistic) %>%
        select(- result)
    } else {
      stop("please supply a valid test")
    }

    taxa_mean_rel_abundances <-
      mean_rel_abundances %>%
      mutate_at(condition, ~ str_c("mean_rel_abundance_in", ., sep = "_")) %>%
      spread(value = mean_rel_abundance, key = condition) %>%
      left_join(taxa_test, by = "taxon_id")

  } else {

    taxa_mean_rel_abundances <-
      mean_rel_abundances %>%
      mutate_at(condition, ~ str_c("mean_rel_abundance_in", ., sep = "_")) %>%
      spread(value = mean_rel_abundance, key = condition)

  }

  # cleanup
  if (exists("rel_abundance_tmp")) ta$counts$rel_abundance <- NULL

  ta %>%
    purrr::modify_at("taxa", left_join, taxa_mean_rel_abundances, by = "taxon_id")

}
