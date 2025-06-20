#' (Re)classify amplicon sequences
#'
#' This function requires the DADA2 package to be installed.
#'
#' `classify_taxa()` will (re)classify either all or a subset of the taxa, given
#' that a variable is present in the taxon table that contains (representative)
#' sequences of the taxa.
#'
#' Ranks can be supplied as a named integer vector, where the names represent
#' taxonomic ranks and the integers represent positions of these ranks in the
#' taxonomy strings present in the reference database. Ranks can also be
#' supplied as just a character vector with the rank names; in that case, it is
#' assumed that the database taxonomy string follows the default order
#' (domain, phylum, class, order, family, genus, species).
#' If no ranks are supplied, taxa will be (re)classified at all default ranks.
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
#' @examples
#'
#' # we create a mock database
#' x <- c(
#' ">Level1;Level2;Level3;Level4;Level5;Level6;",
#' "ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCATAGGATAC", 
#' ">Level1;Level2;Level3;Level4;Level5;",
#' "CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC"
#' )
#' # and write it to a file
#' write(x, file="tmp-db.fna")
#'
#' urt_reclass <- urt %>%
#' # filter out samples to save time
#' filter_samples(sample_id %in% c("s1","s2")) %>%
#' classify_taxa(
#'   "tmp-db.fna", n_ranks = 6, # the mock database is used here
#'   ranks=c("kingdom","phylum", "class", "order", "family", "genus")
#' )
#' # remove the temp file
#' unlink("tmp-db.fna")
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

#' Add sensible taxon name to taxon table
#'
#' `add_taxon_name()` creates sensible taxa names by
#' -under default conditions- combining the genus name with a number.
#' The number is only added if there is more than one taxon of that genus.
#' The number indicates the rank of abundance, with 1 indicating
#' the taxon with the highest mean relative abundance within the genus.
#' If genus classification is not available
#' the next most detailed taxonomic rank which is available is used.
#' The sensible taxon name is added to the taxon table
#' under the column name `taxon_name`.
#'
#'
#' @param ta A tidytacos object.
#' @param method The method by which to arrange the taxon names. Currently only
#'   mean_rel_abundance.
#' @param include_species Whether to include the species name or not.
#' @return A tidytacos object.
#' @examples 
#' urt_g <- urt %>% add_taxon_name()
#' # add the species name if present (which is often uncertain in amplicon data)
#' urt_s <- urt %>% add_taxon_name(include_species = TRUE)
#' @importFrom stats na.omit
#' @family taxa-modifiers
#' @export
add_taxon_name <- function(
  ta, method = "mean_rel_abundance", include_species = FALSE
) {

  mean_rel_abundance <- NULL
  best_classification <- arrange_by_me <- n_taxa <- taxon_number <- NULL

  if (method == "mean_rel_abundance") {

    # if mean_rel_abundance not present: add temporarily
    mean_rel_ab_tmp <- ! "mean_rel_abundance" %in% names(ta$taxa)
    if (mean_rel_ab_tmp) ta <- add_mean_rel_abundance(ta)

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
            classification <- as.character(list(...))
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
    mutate(taxon_number = if_else(
      n_taxa > 1, as.character(seq_len(n())), ""
    )) %>%
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

#' Add taxon color for visualization.
#'
#' `add_rel_abundance()` determines the most abundant taxa and assigns
#' them a color for consistent color codes of each taxon in visualizations.
#' A rank can be supplied to aggregate colors higher than the current rank.
#'
#' @param ta A tidytacos object.
#' @param method The method by which to arrange the taxon names. Currently only
#'   mean_rel_abundance or dominance.
#' @param n An integer denoting the amount of most abundant taxa to display.
#'   Capacity at 12.
#' @param samples An optional vector of sample_id's of interest.
#' @param taxa An optional vector of taxon_id's of interest.
#' @param rank An optional rank to aggregate taxa on.
#' @param threshold_dominance An optional threshold for the dominance method.
#' @return A tidytacos object.
#' @family taxa-modifiers
#' @examples
#' # display the 5 most abundant taxa at genus lvl
#' urt %>% add_taxon_name_color(n=5, rank='genus') %>% tacoplot_stack()
#' @export
add_taxon_name_color <- function(
  ta, method = "mean_rel_abundance", n = 12, samples = NULL, taxa = NULL,
  rank = NULL, threshold_dominance = NULL
  ) {

    arrange_by_me <- NULL
    valid_methods <- c("mean_rel_abundance", "dominance")

  # aggregate rank if asked for
  if(! is.null(rank)) {
    ta <- aggregate_taxa(ta, rank = rank)
  }

  # if taxon_name not present: add temporarily
  taxon_name_tmp <- ! "taxon_name" %in% names(ta$taxa)
  if (taxon_name_tmp) ta <- add_taxon_name(ta)

  if (method %in% valid_methods) {

    # if mean_rel_abundance not present: add temporarily
    mean_rel_ab_tmp <- ! "mean_rel_abundance" %in% names(ta$taxa)
    if (mean_rel_ab_tmp) ta <- add_mean_rel_abundance(ta)

    ta <- mutate_taxa(ta, arrange_by_me = mean_rel_abundance)

  } else {
    # throw error if method unknown
    stop("method unknown, please choose from: ",
         paste(valid_methods, collapse = ", "))
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

  if (! is.null(threshold_dominance)) {
    warning("threshold_dominance is set, using dominance method")
    method <- "dominance"
  } else {
    threshold_dominance <- .5
  }

  if (method == "dominance") {
    # if dominant taxa not present: add temporarily
    if (!"dominant_taxon" %in% names(ta_subset$samples)) {
      ta_subset <- add_dominant_taxa(ta_subset,
        threshold_dominance = threshold_dominance)
    } 
    dom_taxa <- ta_subset$samples$dominant_taxon %>%
      na.omit() %>%
      unique()

    for (col in c("taxon_id", ta_subset %>% rank_names())) {
      if (all(dom_taxa %in% ta_subset$taxa[[col]])) {
        ta_subset <- filter_taxa(ta_subset, .data[[col]] %in% dom_taxa)
        break
      }
    }
      
  }

  # extract taxon names to visualize, in order
  levels <-
    ta_subset$taxa %>%
    arrange(desc(arrange_by_me)) %>%
    pull(taxon_name) %>%
    `[`(1:(n-1)) %>%
    sort()
  levels <- append(levels, "Other taxa", after=0)

  # add taxon_name_color factor to taxa table
  ta$taxa <-
    ta$taxa %>%
    mutate(taxon_name_color = if_else(taxon_name %in% levels, taxon_name, "Other taxa")) %>%
    mutate(taxon_name_color = factor(taxon_name_color, levels = levels))

  # cleanup
  if (taxon_name_tmp) ta$taxa$taxon_name <- NULL
  if (exists("mean_rel_ab_tmp")) ta$taxa$mean_rel_abundance <- NULL
  ta$taxa$arrange_by_me <- NULL

  # return ta object
  ta

}

#' Apply the taxon QC method of Jervis-Bardy
#'
#' `add_jervis_bardy()` calculates the spearman correlation between
#' relative abundance and sample DNA concentration,
#' for each taxon and adds the correlation metric and p-value
#' to the taxa table under the column names "jb_cor" and "jb_p", respectively.
#' If taxa show a distribution that is negatively correlated with
#' DNA concentration, it indicates their potential as contaminants.
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
#' @examples
#' library(dplyr)
#' # filter out blank samples
#' plants <- leaf %>%
#'   filter_samples(Plant != "Blank")
#' # assume Leafweight is a proxy for DNA concentration of the sample
#' plants_jb <- plants %>%
#'   add_jervis_bardy(dna_conc = Leafweight)
#'
#' # we can do this in one step!
#' plants_jb <- leaf %>%
#'   add_jervis_bardy(
#'     dna_conc = Leafweight,
#'     sample_condition = Plant != "Blank"
#')
#'
#' # show the negative correlations
#' plants_jb$taxa %>%
#'   select(taxon_id, starts_with("jb_")) %>%
#'   filter(jb_cor < 0) %>%
#'   arrange(jb_p)
#' @return A tidytacos object with the Jervis-Bardy metrics
#' added to the taxa table.
#' @export
add_jervis_bardy <- function(ta, dna_conc,
  sample_condition = TRUE, min_pres = 3
) {

  jb <- . <- NULL
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

#' Add taxon prevalences to the taxon table
#'
#' `add_prevalence()` calculates taxon prevalences
#' (overall or per condition) and adds it to the taxa table
#' under the column name "prevalence".
#' Prevalence can be expressed as the number of samples where a taxon occurs
#' or the ratio of samples where a taxon occurs and the total amount of samples.
#'
#' If 'condition' is specified, the prevalences will
#' be calculated separately for each group defined by the condition variable.
#' This variable should be present in the sample table.
#'
#' If `condition` is specified, differential prevalence testing
#' can be performed by setting the `fisher_test` argument.
#' Options are F (default) or T. When set to T, significance of
#' differential prevalence will be added to the taxa table
#' under column name `fisher_p`.
#'
#' Condition should be a categorical variable present in the samples table.
#' Supply condition as a string.
#' @importFrom stats fisher.test
#' @param ta A tidytacos object.
#' @param condition A categorical variable (string).
#' @param relative Whether to use relative occurrences.
#' @param fisher_test Whether to perform a fisher test and
#' add the p-values of the test to the taxa table.
#' @return A tidytacos object.
#' @examples
#' # add prevalences of all taxa
#' urtp <- urt %>% add_prevalence()
#' urtp$taxa %>% dplyr::select(taxon_id, prevalence)
#'
#' # add prevalences and fisher test for location
#' urtpf <- urt %>%
#'   add_prevalence(condition="location", fisher_test=TRUE, relative=TRUE)
#' urtpf$taxa %>%
#'   dplyr::select(taxon_id, prevalence_in_N, prevalence_in_NF, fisher_p)
#' @family taxa-modifiers
#' @export
add_prevalence <- function(
  ta, condition = NULL, relative = FALSE, fisher_test = FALSE
) {

  # remove any pre-existing result columns
  ta$taxa <- ta$taxa %>%
    select(- starts_with("prevalence"), - starts_with("fisher"))

  fisher <- prevalence <- . <- NULL
  prev <- "prevalence"
  prev_in <- "prevalence_in"
  if (is.null(condition)) {

    taxa_prevalences <-
      prevalences(ta, condition = condition)

  } else if (fisher_test) {

    prevalences <-
      prevalences(ta, condition = condition, pres_abs = TRUE)

    condition_sym <- sym(condition)

    taxa_fisher <-
      prevalences %>%
      group_by(taxon_id) %>%
      arrange(!! condition_sym, presence) %>%
      do(
        fisher = c(.$n) %>%
          matrix(ncol = 2, byrow = TRUE) %>%
          fisher.test()
      ) %>%
      mutate(fisher_p = fisher$p.value) %>%
      select(- fisher)

    taxa_prevalences <-
      prevalences %>%
      filter(presence == "present") %>%
      select(taxon_id, !! condition_sym, prevalence = n) %>%
      mutate_at(condition, ~ str_c(prev_in, ., sep = "_")) %>%
      spread(value = prev, key = condition) %>%
      left_join(taxa_fisher, by = "taxon_id")

  } else {

    taxa_prevalences <-
      prevalences(ta, condition = condition) %>%
      mutate_at(condition, ~ str_c(prev_in, ., sep = "_")) %>%
      spread(value = prev, key = condition)

  }

  if (relative && is.null(condition)) {

    taxa_prevalences <-
      taxa_prevalences %>%
      mutate(prevalence = prevalence / nrow(ta$samples))

  }

  if (relative && ! is.null(condition)) {

    condition_sym <- sym(condition)

    conditions <-
      ta$samples %>%
      count(!! condition_sym)

    for (con_ix in seq_len(nrow(conditions))) {

      con <- conditions[[condition]][con_ix]
      n_samples <- conditions[["n"]][con_ix]
      taxa_prevalences[[str_c(prev_in, con, sep="_")]] <-
        taxa_prevalences[[str_c(prev_in, con, sep="_")]] / n_samples

    }

  }

  ta %>%
    purrr::modify_at("taxa", left_join, taxa_prevalences, by = "taxon_id")

}

#' Add average relative abundances to taxa table
#'
#' `add_mean_rel_abundance()` adds mean relative abundance values
#' for each taxon to the taxa table, overall or per sample group.
#'
#' If `condition` is specified, the mean relative abundances will be calculated
#' separately for each group defined by the condition variable. This variable
#' should be present in the sample table.
#'
#' If `condition` is specified, differential abundance testing can be performed
#' by setting the `test` argument. Options are `NULL` (default), `"wilcox"` or
#' `"t-test"`.
#'
#' @importFrom stats t.test wilcox.test
#' @param ta A tidytacos object.
#' @param condition A condition variable (character).
#' @param test Differential abundance test to perform.
#'
#' @return A tidytacos object
#' @family taxa-modifiers
#'
#' @export
add_mean_rel_abundance <- function(ta, condition = NULL, test = NULL) {

  result <- mean_rel_abundance <- . <- NULL
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
    purrr::modify_at(
      "taxa", left_join, taxa_mean_rel_abundances, by = "taxon_id"
    )

}
